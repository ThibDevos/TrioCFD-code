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

#include <vector>

#include <Postraitement_DT.h>
#include <SFichier.h>

Implemente_instanciable(Postraitement_DT, "Postraitement_DT", Postraitement);

/**
 * @brief Reads the input stream and initializes the output file.
 *
 * This method reads the necessary parameters, removes the file extension
 * from `nom_fich_`, and writes the header information to the output file.
 *
 * @param is Input stream.
 * @return Entree& The modified input stream.
 */
Entree& Postraitement_DT::readOn(Entree& is) {
	// read input has standard Postraitement
	Postraitement::readOn(is);

	// Remove file extension from nom_fich_ that was automatically added
	std::string file_name(nom_fich_);
	size_t last_dot_pos = file_name.find_last_of('.');
	if (last_dot_pos != std::string::npos)
		nom_fich_ = file_name.substr(0, last_dot_pos);

	// get problem
	const Probleme_base& problem = ref_cast(Probleme_base, mon_probleme.valeur());

	// Open file
	SFichier file;
	file.ouvrir(nom_fich(), ios::out);

	// allocate memory for a text buffer
	size_t buffer_size = static_cast<size_t>(std::max(100, column_width + column_gap + 1));
	char* buffer = new char[buffer_size];

	// Write description
	file << "File generated from Postraitement_DT class.\n\n";

	file << format("Description:\n", "\033[1m");
	file << "This file contains the time steps computed for each equations, has well has the minimum, maximum and effective time steps.\n";
	file << "All values are in seconds.\n\n";

	// Write equations
	file << format("Equations [type and name]:\n", "\033[1m");
	for (int index = 0; index < problem.nombre_d_equations(); index++) {
		const Nom& equation_type = problem.equation(index).que_suis_je();
		const Nom& equation_name = problem.equation(index).le_nom();

		// Print equation type and name using tab
		snprintf(buffer , buffer_size, "\t%-47s\t%-15s\n",
			&(*equation_type),
			&(*equation_name));
		file << buffer;
	}
	file << "\n";

	// print the color code used
	if (formatting_flag) {
		file << "\033[1mColor code:\n";
		
		// get formatters
		std::string formatter_minimum = get_formatter(Postraitement_DT::Extrema::Min);
		std::string formatter_maximum = get_formatter(Postraitement_DT::Extrema::Max);
		std::string formatter_effective = get_formatter(Postraitement_DT::Extrema::Effective);
		std::string formatter_reset = "\033[0m";

		// write
		file << formatter_minimum << "\tMinimum" << formatter_reset
		     << ":   The minimum value among all time steps for the given time.\n"
		     << formatter_maximum << "\tMaximum" << formatter_reset
			 << ":   The maximum value among all time steps for the given time.\n"
			 << formatter_effective << "\tEffective" << formatter_reset
			 << ": The time step value actually used for computation at the given time."
			 << "\033[0m\n\n";
	}

	// write the header
	if (formatting_flag)
		file << "\033[1m"; // bold

	write_text_column(buffer, buffer_size, "Time");
	file << buffer;

	for (int index = 0; index < mon_probleme.valeur().nombre_d_equations(); index++) {
		const Nom& equation_name = problem.equation(index).le_nom();
		write_text_column(buffer, buffer_size, &(*equation_name));
		file << buffer;
	}

	write_text_column(buffer, buffer_size, "Minimum");
	file << buffer;

	write_text_column(buffer, buffer_size, "Maximum");
	file << buffer;

	write_text_column(buffer, buffer_size, "Effective");
	file << buffer;

	if (formatting_flag)
		file << "\033[0m"; // normal

	file << "\n";

	// free the buffer 
	delete[] buffer;

	// close the file
	file.close();

	return is;
}

/**
 * @brief Prints the object to the output stream.
 *
 * This method currently does nothing and returns the output stream as is.
 *
 * @param os Output stream.
 * @return Sortie& The unchanged output stream.
 */
Sortie& Postraitement_DT::printOn(Sortie& os) const {
	// do nothing
	return os;
}

/**
 * @brief Sets the parameters of the object.
 *
 * This method delegates parameter setting to the base class `Postraitement`.
 *
 * @param param Reference to the parameter object.
 */
void Postraitement_DT::set_param(Param& param) {
  Postraitement::set_param(param);

  param.ajouter("column_width",       &column_width);
  param.ajouter("column_gap",         &column_gap);
  param.ajouter("number_of_decimals", &number_of_decimals);
  param.ajouter("formatting",         &formatting_flag);
}

/**
 * @brief Appends a new line to the output file with the estimated time steps.
 * 
 * This method writes the current simulation time followed by the time step values
 * for each equation of the problem to the output file.
 */
void Postraitement_DT::postraiter(int) {
	
	// get problem and time scheme
	const Probleme_base& problem = ref_cast(Probleme_base, mon_probleme.valeur());
	const Schema_Temps_base& time_scheme = ref_cast(Schema_Temps_base, problem.schema_temps());

	// Open file
	SFichier file;
	file.ouvrir(nom_fich(), ios::app);

	// allocate memory for a text buffer
	size_t buffer_size = column_width + column_gap + 1;
	char* buffer = new char[buffer_size];

	// Get current simulation time and write it to the file 
	double current_time = time_scheme.temps_courant();
	write_float_column(buffer, buffer_size, current_time);
	file << buffer;

	// get all time steps (equations, minimum, maximum and effective time steps)
	std::vector<double> time_steps;

	for (int index = 0; index < problem.nombre_d_equations(); index++) {
		double dt = problem.equation(index).calculer_pas_de_temps();
		time_steps.push_back(dt);
	}

	time_steps.push_back(time_scheme.pas_temps_min());
	time_steps.push_back(time_scheme.pas_temps_max());

	double effective_time_step = problem.calculer_pas_de_temps();
	time_steps.push_back(effective_time_step);

	// find the minimum and maximum time step values
	auto [min_it, max_it] = std::minmax_element(time_steps.begin(), time_steps.end());

	// Write all time steps
	for (double time_step: time_steps) {
		// format differently the output if its a minimum or maximum value
		std::string formatter;
		if       (time_step == effective_time_step)
			formatter = get_formatter(Postraitement_DT::Extrema::Effective);
		else if (time_step == *min_it)
			formatter = get_formatter(Postraitement_DT::Extrema::Min);
		else if (time_step == *max_it)
			formatter = get_formatter(Postraitement_DT::Extrema::Max);

		// write the value in the buffer
		write_float_column(buffer, buffer_size, time_step);

		// write to the files with formatter
		file << format(buffer, formatter);
	}
	file << "\n";

	// free the buffer 
	delete[] buffer;

	// close the file
	file.close();
}

/**
 * @brief Helper function to get the formatting to apply to a text output.
 *
 * @param extrema Extrema value that condition the formating.
 * @return std::string A prefix that can be use has formatting string.
 */
std::string Postraitement_DT::get_formatter(const Postraitement_DT::Extrema& extrema) const {
	if      (extrema == Postraitement_DT::Extrema::Min)
		return "\033[31;1m"; // red bold
	else if (extrema == Postraitement_DT::Extrema::Max)
		return "\033[32;1m"; // red bold
	else if (extrema == Postraitement_DT::Extrema::Effective)
		return "\033[33;1m"; // red bold
	else
		return ""; // normal
}

/**
 * @brief apply or not the formatter to the string depending on the formatting_flag.
 *
 * @param input Input string.
 * @param formatter Formatter string.
 * @return std::string Formatted string.
 */
std::string Postraitement_DT::format(
		const std::string input, const std::string formatter) const {
	
	if (formatting_flag)
		return formatter + input + "\033[0m";
	else
		return input;
}

/**
 * @brief Writes a floating point value to a character buffer.
 * 
 * This function formats a floating point value into scientific notation if 
 * greater than 1, ensuring a fixed column width with a trailing gap.
 * 
 * @tparam Float Floating point type (e.g., float, double).
 * @param buffer Character buffer to store the formatted output.
 * @param size Size of the buffer.
 * @param value Floating point value to format.
 */
template <class Float>
void Postraitement_DT::write_float_column(
		char buffer[], size_t size, const Float value) const {

	// check template type
	static_assert(std::is_floating_point<Float>::value, "write_float_column requires a floating-point type (float or double).");

	// format the value using snprintf and add a gap (the column gap) at the end
	snprintf(buffer, size, "%-*.*e%*s",
		column_width, number_of_decimals, value,
		column_gap, "");
}

/**
 * @brief Writes a text value to a character buffer.
 * 
 * This function formats a text string into a fixed-width column with a trailing gap.
 * 
 * @param buffer Character buffer to store the formatted output.
 * @param size Size of the buffer.
 * @param text Text string to format.
 */
void Postraitement_DT::write_text_column(
		char buffer[], size_t size, const char text[]) const {

	// format the text using snprintf and add a gap (the column gap) at the end
	snprintf(buffer, size, "%-*s%*s",
		column_width, text,
		column_gap, "");
}
