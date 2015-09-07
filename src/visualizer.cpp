/***************************************************************************
 *   Copyright (C) 2008-2014 by Andrzej Rybczak                            *
 *   electricityispower@gmail.com                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.              *
 ***************************************************************************/

#include "visualizer.h"

#ifdef ENABLE_VISUALIZER

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/math/constants/constants.hpp>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <fcntl.h>

#include "global.h"
#include "settings.h"
#include "status.h"
#include "statusbar.h"
#include "title.h"
#include "screen_switcher.h"
#include "status.h"
#include "enums.h"

using Global::MainStartY;
using Global::MainHeight;

// Visualizer object
Visualizer *myVisualizer;

namespace {

// FPS ???
const int fps = 25;

// toColor: a scaling function for coloring. For numbers 0 to max this function returns
// a coloring from the lowest color to the highest, and colors will not loop from 0 to max.
const NC::Color &toColor(size_t number, size_t max, bool wrap = true)
{
	const auto colors_size = Config.visualizer_colors.size();
	const auto index = (number * colors_size) / max;
	return Config.visualizer_colors[
		wrap ? index % colors_size : std::min(index, colors_size-1)
	];
}

}

Visualizer::Visualizer()
: Screen(NC::Window(0, MainStartY, COLS, MainHeight, "", NC::Color::Default, NC::Border()))
{
	ResetFD();
	// this is the number of samples it gets from MPD (every second?)
	m_samples = 44100/fps;
	// with stereo, there's double the samples
	if (Config.visualizer_in_stereo)
		m_samples *= 2;
	
	// I don't know wtf this all is
#	ifdef HAVE_FFTW3_H
	m_fftw_results = m_samples/2+1;
	m_freq_magnitudes.resize(m_fftw_results);
	m_fftw_input = static_cast<double *>(fftw_malloc(sizeof(double)*m_samples));
	m_fftw_output = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*m_fftw_results));
	m_fftw_plan = fftw_plan_dft_r2c_1d(m_samples, m_fftw_input, m_fftw_output, FFTW_ESTIMATE);
#	endif // HAVE_FFTW3_H
}

// when switching to the visualizer within ncmpcpp
void Visualizer::switchTo()
{
	// call the constructor
	SwitchTo::execute(this);
	// clear the window
	w.clear();
	// TODO: need to figure out what this does
	SetFD();
	// TODO: need to figure out what this is
	m_timer = boost::posix_time::from_time_t(0);
	drawHeader();
}

// when the ncurses window is resized
void Visualizer::resize()
{
	size_t x_offset, width;
	// get the new x and y dimensions
	getWindowResizeParams(x_offset, width);
	// resize the ncurses window
	w.resize(width, MainHeight);
	// move the stuff in the window to fit the new window size ???
	w.moveTo(x_offset, MainStartY);
	// reset the resize flag
	hasToBeResized = 0;
}

// getter for the visualizer pane title
std::wstring Visualizer::title()
{
	// TODO: what does the L mean
	return L"Music visualizer";
}

void Visualizer::update()
{
	// if fifo doesn't exist
	if (m_fifo < 0)
		return;

	// PCM in format 44100:16:1 (for mono visualization) and
	// 44100:16:2 (for stereo visualization) is supported.
	int16_t buf[m_samples]; // create array
	// read from the fifo buffer
	ssize_t data = read(m_fifo, buf, sizeof(buf));
	if (data < 0) // no data available in fifo
		return;

	// TODO: if the visualizer is out of sync ???
	if (m_output_id != -1 && Global::Timer - m_timer > Config.visualizer_sync_interval)
	{
		Mpd.DisableOutput(m_output_id);
		usleep(50000); // TODO: sleep for 50 seconds???
		Mpd.EnableOutput(m_output_id);
		m_timer = Global::Timer; // put the timer in sync
	}

	void (Visualizer::*draw)(int16_t *, ssize_t, size_t, size_t);
	void (Visualizer::*drawStereo)(int16_t *, int16_t *, ssize_t, size_t);
#	ifdef HAVE_FFTW3_H // TODO: fast fourier transform... what is W3_H???
	if (Config.visualizer_type == VisualizerType::Spectrum)
	{
		draw = &Visualizer::DrawFrequencySpectrum;
		drawStereo = &Visualizer::DrawFrequencySpectrumStereo;
	}
	else
#	endif // HAVE_FFTW3_H
	if (Config.visualizer_type == VisualizerType::WaveFilled)
	{
		draw = &Visualizer::DrawSoundWaveFill;
		drawStereo = &Visualizer::DrawSoundWaveFillStereo;
	}
	else if (Config.visualizer_type == VisualizerType::Ellipse) // TODO: how do you enable this???
	{
		draw = &Visualizer::DrawSoundEllipse;
		drawStereo = &Visualizer::DrawSoundEllipseStereo;
	}
	else
	{
		draw = &Visualizer::DrawSoundWave;
		drawStereo = &Visualizer::DrawSoundWaveStereo;
	}

	// number of samples is the size of the data/16
	// therefore each sample is 16 bits
	const ssize_t samples_read = data/sizeof(int16_t);
	if (Config.visualizer_sample_multiplier == 1.0)
	{
		m_auto_scale_multiplier += 1.0/fps;
		std::for_each(buf, buf+samples_read, [this](int16_t &sample) {
			double scale = std::numeric_limits<int16_t>::min();
			scale /= sample;
			scale = fabs(scale); // TODO: what is this? Absolute value???
			if (scale < m_auto_scale_multiplier)
				m_auto_scale_multiplier = scale;
		});
	}
	// for each sample in the buffer, scale the sample
	std::for_each(buf, buf+samples_read, [this](int16_t &sample) {
		int32_t tmp = sample;
		if (Config.visualizer_sample_multiplier != 1.0)
			tmp *= Config.visualizer_sample_multiplier;
		else if (m_auto_scale_multiplier <= 50.0) // limit the auto scale
			tmp *= m_auto_scale_multiplier;
		if (tmp < std::numeric_limits<int16_t>::min())
			sample = std::numeric_limits<int16_t>::min();
		else if (tmp > std::numeric_limits<int16_t>::max())
			sample = std::numeric_limits<int16_t>::max();
		else
			sample = tmp;
	});

	// clear the screen
	w.clear();
	if (Config.visualizer_in_stereo)
	{
		// TODO: what does auto mean???
		// the number of samples is halfed because half are for the left channel and half are for the right
		auto chan_samples = samples_read/2;
		// makes arrays for left and right audio channels
		int16_t buf_left[chan_samples], buf_right[chan_samples];
		// alternate samples between the left buffer and right buffer
		for (ssize_t i = 0, j = 0; i < samples_read; i += 2, ++j)
		{
			buf_left[j] = buf[i];
			buf_right[j] = buf[i+1];
		}
		size_t half_height = w.getHeight()/2;

		// draw the visualizer
		(this->*drawStereo)(buf_left, buf_right, chan_samples, half_height);
	}
	else
	{
		// draw the visualizer
		(this->*draw)(buf, samples_read, 0, w.getHeight());
	}
	// refresh the window
	w.refresh();
}

// TODO: I don't know what this does or when it does it, but I know why it does
int Visualizer::windowTimeout()
{
	if (m_fifo >= 0 && Status::State::player() == MPD::psPlay)
		return 1000/fps;
	else
		return Screen<WindowType>::windowTimeout();
}

// TODO: what does this line separate???
/**********************************************************************/

// sound wave visualizer
void Visualizer::DrawSoundWave(int16_t *buf, ssize_t samples, size_t y_offset, size_t height)
{
	// initialize variables
	const size_t half_height = height/2;
	const size_t base_y = y_offset+half_height;
	const size_t win_width = w.getWidth();
	const int samples_per_column = samples/win_width;

	// too little samples
	if (samples_per_column == 0)
		return;

	// TODO: what is auto???
	// TODO: what does all the <<s do??
	auto draw_point = [&](size_t x, int32_t y) {
		w << NC::XY(x, base_y+y)
		<< toColor(std::abs(y), half_height, false)
		<< Config.visualizer_chars[0]
		<< NC::Color::End;
	};

	int32_t point_y, prev_point_y = 0;
	// for each x value, calculate the y
	for (size_t x = 0; x < win_width; ++x)
	{
		point_y = 0;
		// calculate mean from the relevant points
		for (int j = 0; j < samples_per_column; ++j)
			point_y += buf[x*samples_per_column+j];
		point_y /= samples_per_column;
		// normalize it to fit the screen
		point_y *= height / 65536.0;

		// draw the point on the screen
		draw_point(x, point_y);

		// if the gap between two consecutive points is too big,
		// intermediate values are needed for the wave to be watchable.
		if (x > 0 && std::abs(prev_point_y-point_y) > 1)
		{
			const int32_t half = (prev_point_y+point_y)/2;
			if (prev_point_y < point_y)
			{
				for (auto y = prev_point_y; y < point_y; ++y)
					draw_point(x-(y < half), y);
			}
			else
			{
				for (auto y = prev_point_y; y > point_y; --y)
					draw_point(x-(y > half), y);
			}
		}
		prev_point_y = point_y;
	}
}

// if it's in stereo, just draw two sound waves
void Visualizer::DrawSoundWaveStereo(int16_t *buf_left, int16_t *buf_right, ssize_t samples, size_t height)
{
	DrawSoundWave(buf_left, samples, 0, height);
	DrawSoundWave(buf_right, samples, height, w.getHeight() - height);
}

/**********************************************************************/

// DrawSoundWaveFill: This visualizer is very similar to DrawSoundWave, but instead of
// a single line the entire height is filled. In stereo mode, the top half of the screen
// is dedicated to the right channel, the bottom the left channel.
void Visualizer::DrawSoundWaveFill(int16_t *buf, ssize_t samples, size_t y_offset, size_t height)
{
	// if right channel is drawn, bars descend from the top to the bottom
	const bool flipped = y_offset > 0;
	const size_t win_width = w.getWidth();
	const int samples_per_column = samples/win_width;

	// too little samples
	if (samples_per_column == 0)
		return;

	int32_t point_y;
	for (size_t x = 0; x < win_width; ++x)
	{
		point_y = 0;
		// calculate mean from the relevant points
		for (int j = 0; j < samples_per_column; ++j)
			point_y += buf[x*samples_per_column+j];
		point_y /= samples_per_column;
		// normalize it to fit the screen
		point_y = std::abs(point_y);
		point_y *= height / 32768.0;

		for (int32_t j = 0; j < point_y; ++j)
		{
			size_t y = flipped ? y_offset+j : y_offset+height-j-1;
			// TODO: what do all the <<s do???
			w << NC::XY(x, y)
			<< toColor(j, height)
			<< Config.visualizer_chars[1]
			<< NC::Color::End;
		}
	}
}

// just draw two soundwaves if the song is in stereo
void Visualizer::DrawSoundWaveFillStereo(int16_t *buf_left, int16_t *buf_right, ssize_t samples, size_t height)
{
	DrawSoundWaveFill(buf_left, samples, 0, height);
	DrawSoundWaveFill(buf_right, samples, height, w.getHeight() - height);
}

/**********************************************************************/

// can't figure out how to activate this
// draws the sound wave as an ellipse with origin in the center of the screen
void Visualizer::DrawSoundEllipse(int16_t *buf, ssize_t samples, size_t, size_t height)
{
	// initialize variables
	// these are probably used to calculate where the center of the circle should be
	const size_t half_width = w.getWidth()/2;
	const size_t half_height = height/2;

	// make it so that the loop goes around the ellipse exactly once
	// TODO: what is boost???
	const double deg_multiplier = 2*boost::math::constants::pi<double>()/samples;

	int32_t x, y;
	double radius, max_radius;
	// for each of the samples
	for (ssize_t i = 0; i < samples; ++i)
	{
		x = half_width * std::cos(i*deg_multiplier);
		y = half_height * std::sin(i*deg_multiplier);
		max_radius = sqrt(x*x + y*y);

		// calculate the distance of the sample from the center,
		// where 0 is the center of the ellipse and 1 is its border
		radius = std::abs(buf[i]);
		radius /= 32768.0;

		// appropriately scale the position
		x *= radius;
		y *= radius;

		w << NC::XY(half_width + x, half_height + y)
		<< toColor(sqrt(x*x + y*y), max_radius, false)
		<< Config.visualizer_chars[0]
		<< NC::Color::End;
	}
}

// DrawSoundEllipseStereo: This visualizer only works in stereo. The colors form concentric
// rings originating from the center (width/2, height/2). For any given point, the width is
// scaled with the left channel and height is scaled with the right channel. For example,
// if a song is entirely in the right channel, then it would just be a vertical line.
//
// Since every font/terminal is different, the visualizer is never a perfect circle. This
// visualizer assume the font height is twice the length of the font's width. If the font
// is skinner or wider than this, instead of a circle it will be an ellipse.
void Visualizer::DrawSoundEllipseStereo(int16_t *buf_left, int16_t *buf_right, ssize_t samples, size_t half_height)
{
	const size_t width = w.getWidth();
	const size_t left_half_width = width/2;
	const size_t right_half_width = width - left_half_width;
	const size_t top_half_height = half_height;
	const size_t bottom_half_height = w.getHeight() - half_height;

	// Makes the radius of each ring be approximately 2 cells wide.
	const int32_t radius = 2*Config.visualizer_colors.size();
	int32_t x, y;
	for (ssize_t i = 0; i < samples; ++i)
	{
		x = buf_left[i]/32768.0 * (buf_left[i] < 0 ? left_half_width : right_half_width);
		y = buf_right[i]/32768.0 * (buf_right[i] < 0 ? top_half_height : bottom_half_height);

		// The arguments to the toColor function roughly follow a circle equation where
		// the center is not centered around (0,0). For example (x - w)^2 + (y-h)+2 = r^2
		// centers the circle around the point (w,h). Because fonts are not all the same
		// size, this will not always generate a perfect circle.
		w << toColor(sqrt(x*x + 4*y*y), radius)
		<< NC::XY(left_half_width + x, top_half_height + y)
		<< Config.visualizer_chars[1]
		<< NC::Color::End;
	}
}

/**********************************************************************/

#ifdef HAVE_FFTW3_H
void Visualizer::DrawFrequencySpectrum(int16_t *buf, ssize_t samples, size_t y_offset, size_t height)
{
	// if right channel is drawn, bars descend from the top to the bottom
	const bool flipped = y_offset > 0;

	// copy samples to fftw input array
	for (unsigned i = 0; i < m_samples; ++i)
		m_fftw_input[i] = i < samples ? buf[i] : 0;
	fftw_execute(m_fftw_plan);

	// count magnitude of each frequency and scale it to fit the screen
	for (size_t i = 0; i < m_fftw_results; ++i)
		m_freq_magnitudes[i] = sqrt(
			m_fftw_output[i][0]*m_fftw_output[i][0]
		+	m_fftw_output[i][1]*m_fftw_output[i][1]
		)/2e4*height;

	const size_t win_width = w.getWidth();
	// cut bandwidth a little to achieve better look
	const double bins_per_bar = m_fftw_results/win_width * 7/10;
	double bar_height;
	size_t bar_bound_height;
	for (size_t x = 0; x < win_width; ++x)
	{
		bar_height = 0;
		for (int j = 0; j < bins_per_bar; ++j)
			bar_height += m_freq_magnitudes[x*bins_per_bar+j];
		// buff higher frequencies
		bar_height *= log2(2 + x) * 100.0/win_width;
		// moderately normalize the heights
		bar_height = pow(bar_height, 0.5);

		bar_bound_height = std::min(std::size_t(bar_height/bins_per_bar), height);
		for (size_t j = 0; j < bar_bound_height; ++j)
		{
			size_t y = flipped ? y_offset+j : y_offset+height-j-1;
			w << NC::XY(x, y)
			<< toColor(j, height)
			<< Config.visualizer_chars[1]
			<< NC::Color::End;
		}
	}
}

void Visualizer::DrawFrequencySpectrumStereo(int16_t *buf_left, int16_t *buf_right, ssize_t samples, size_t height)
{
	DrawFrequencySpectrum(buf_left, samples, 0, height);
	DrawFrequencySpectrum(buf_right, samples, height, w.getHeight() - height);
}
#endif // HAVE_FFTW3_H

/**********************************************************************/

void Visualizer::ToggleVisualizationType()
{
	switch (Config.visualizer_type)
	{
		case VisualizerType::Wave:
			Config.visualizer_type = VisualizerType::WaveFilled;
			break;
		case VisualizerType::WaveFilled:
#			ifdef HAVE_FFTW3_H
			Config.visualizer_type = VisualizerType::Spectrum;
#			else
			Config.visualizer_type = VisualizerType::Ellipse;
#			endif // HAVE_FFTW3_H
			break;
#		ifdef HAVE_FFTW3_H
		case VisualizerType::Spectrum:
			Config.visualizer_type = VisualizerType::Ellipse;
			break;
#		endif // HAVE_FFTW3_H
		case VisualizerType::Ellipse:
			Config.visualizer_type = VisualizerType::Wave;
			break;
	}
	Statusbar::printf("Visualization type: %1%", Config.visualizer_type);
}

void Visualizer::SetFD()
{
	if (m_fifo < 0 && (m_fifo = open(Config.visualizer_fifo_path.c_str(), O_RDONLY | O_NONBLOCK)) < 0)
		Statusbar::printf("Couldn't open \"%1%\" for reading PCM data: %2%",
			Config.visualizer_fifo_path, strerror(errno)
		);
}

void Visualizer::ResetFD()
{
	m_fifo = -1;
}

void Visualizer::FindOutputID()
{
	m_output_id = -1;
	if (!Config.visualizer_output_name.empty())
	{
		for (MPD::OutputIterator out = Mpd.GetOutputs(), end; out != end; ++out)
		{
			if (out->name() == Config.visualizer_output_name)
			{
				m_output_id = out->id();
				break;
			}
		}
		if (m_output_id == -1)
			Statusbar::printf("There is no output named \"%s\"", Config.visualizer_output_name);
	}
}

void Visualizer::ResetAutoScaleMultiplier()
{
	m_auto_scale_multiplier = std::numeric_limits<double>::infinity();
}

#endif // ENABLE_VISUALIZER

/* vim: set tabstop=4 softtabstop=4 shiftwidth=4 noexpandtab : */
