/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
#ifndef Postraitement_DT_included
#define Postraitement_DT_included

#include <Postraitement.h>
#include <TRUST_Ref.h>

/*
 * This class when use in the Postraitement block of a .data file allow the creation
 * of an output file giving time-step related informations about the problem being
 * solved.
 * The header of the file is created during readOn and a new line is added everytime
 * the 'postraiter' method is called.
 */

class Postraitement_DT : public Postraitement {

		/////////////////
		// declaration //
		/////////////////

	private:

		Declare_instanciable(Postraitement_DT);

		///////////
		// enums //
		///////////

	public:
		
		// This enum is used in formatting function to specify whether the value
		// that need to be formated needs to be formatted as a the minimum, maximum
		// or as the effective timestep value.
		enum Extrema {
			Min,
			Max,
			Effective,
			None,
		};

		/////////////
		// methods //
		/////////////

		void set_param(Param& param) override;
		void postraiter(int) override;

	private:

		// method that return formatter for output depending of the extrema
		std::string get_formatter(const Extrema& extrema) const;

		// method that apply formatting
		std::string format(const std::string input, const std::string formatter) const;

		// two auxilary methods to format floating points and text values into columns
		template <class Float>
		void write_float_column(
			char buffer[], size_t size, const Float value) const;

		void write_text_column(
			char buffer[], size_t size, const char text[]) const;

		////////////////
		// attributes //
		////////////////

	private:
		
		// text formatting (in number of characters)
		int column_width       = 13;
		int column_gap         = 3;

		// number of decimal number if floating points values
		int number_of_decimals = 7;

		// flag to activate or not formatting
		int formatting_flag    = 1; // by default formatting is active
};

#endif
