/**
 * \file
 *
 * \author Valentin Bruder
 *
 * \copyright Copyright (C) 2018 Valentin Bruder
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <datrawreader.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cassert>

/*
 * DatRawReader::read_files
 */
void DatRawReader::read_files(const std::string dat_file_name)
{
    // check file
    if (!dat_file_name.empty())
        _prop.dat_file_name = dat_file_name;
    else
        throw std::invalid_argument("dat file name must not be empty.");

    try
    {
        this->read_dat(_prop.dat_file_name);      
        this->_raw_data.clear();
        for (const auto &n : _prop.raw_file_names)
            read_raw(n);
    }
    catch (std::runtime_error e)
    {
        throw e;
    }
}


/*
 * DatRawReader::has_data
 */
bool DatRawReader::has_data() const
{
    return !(_raw_data.empty());
}


/*
 * DatRawReader::data
 */
const std::vector<std::vector<char> > & DatRawReader::data() const
{
    if (!has_data())
    {
        throw std::runtime_error("No data available.");
    }
    return _raw_data;
}

/**
 * @brief DatRawReader::properties
 * @return
 */
const Properties &DatRawReader::properties() const
{
    if (!has_data())
    {
        throw std::runtime_error("No properties of volume data set available.");
    }
    return _prop;
}

/**
 * @brief DatRawReader::clearData
 */
void DatRawReader::clearData()
{
    _raw_data.clear();
}


/*
 * DatRawReader::read_dat
 */
void DatRawReader::read_dat(const std::string dat_file_name)
{
    std::ifstream dat_file(dat_file_name);
    std::string line;
    std::vector<std::vector<std::string>> lines;

    // read lines from .dat file and split on whitespace
    if (dat_file.is_open())
    {
        while (std::getline(dat_file, line))
        {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::copy(std::istream_iterator<std::string>(iss),
                      std::istream_iterator<std::string>(),
                      std::back_inserter(tokens));
            lines.push_back(tokens);
        }
        dat_file.close();
    }
    else
    {
        throw std::runtime_error("Could not open .dat file " + dat_file_name);
    }

    for (auto l : lines)
    {
        if (!l.empty())
        {
            std::string name = l.at(0);
            if (name.find("ObjectFileName") != std::string::npos && l.size() > 1u)
            {
                _prop.raw_file_names.clear();
                for (const auto &s : l)
                {
                    if (s.find("ObjectFileName") == std::string::npos)
                        _prop.raw_file_names.push_back(s);
                }
            }
            else if (name.find("Resolution") != std::string::npos && l.size() > 3u)
            {
                for (size_t i = 1; i < l.size(); ++i)
                {
                    _prop.volume_res.at(i - 1) = std::stoul(l.at(i));
                }
            }
            else if (name.find("SliceThickness") != std::string::npos && l.size() > 3u)
            {
                for (size_t i = 1; i < l.size(); ++i)
                {
                    _prop.slice_thickness.at(i - 1) = std::stod(l.at(i));
                }
            }
            if (name.find("Format") != std::string::npos && l.size() > 1u)
            {
                _prop.format = l.at(1);
            }
            if (name.find("Nodes") != std::string::npos && l.size() > 1u)
            {
                _prop.node_file_name = l.at(1);
            }
            if (name.find("TimeSeries") != std::string::npos && l.size() > 1u)
            {
                _prop.time_series = std::stoi(l.at(1));
            }
        }
    }

    // check that values read from the dat file
    if (_prop.raw_file_names.empty())
    {
        throw std::runtime_error("Missing raw file names declaration in " + dat_file_name);
    }
    if (_prop.volume_res.empty())
    {
        throw std::runtime_error("Missing volume resolution declaration in " + dat_file_name);
    }
    if (_prop.slice_thickness.at(0) == -1.0 || _prop.slice_thickness.at(1) == -1.0
            || _prop.slice_thickness.at(2) == -1.0)
    {
        std::cerr << "WARNING: Missing slice thickness declaration in " << dat_file_name
                  << std::endl;
        std::cerr << "Assuming a slice thickness of 1.0 in each dimension." << std::endl;
        _prop.slice_thickness.fill(1.0);
    }
    if (_prop.format.empty())
    {
        std::cerr << "WARNING: Missing format declaration in " << dat_file_name << std::endl;
        std::cerr << "Trying to calculate the format from raw file size and volume resolution."
                  << std::endl;
    }
}

/*
 * DatRawReader::read_raw
 */
void DatRawReader::read_raw(const std::string raw_file_name)
{
    if (raw_file_name.empty())
        throw std::invalid_argument("Raw file name must not be empty.");

    // append .raw file name to .dat file name path
    std::size_t found = _prop.dat_file_name.find_last_of("/\\");
    std::string name_with_path;
    if (found != std::string::npos && _prop.dat_file_name.size() >= found)
    {
        name_with_path = _prop.dat_file_name.substr(0, found + 1) + raw_file_name;
    }
    else
    {
        name_with_path = raw_file_name;
    }

    // use plain old C++ method for file read here that is much faster than iterator
    // based approaches according to:
    // http://insanecoding.blogspot.de/2011/11/how-to-read-in-file-in-c.html
    std::ifstream is(name_with_path, std::ios::in | std::ifstream::binary);
    if (is)
    {
        // get length of file:
        is.seekg(0, is.end);
#ifdef _WIN32
        // HACK: to support files bigger than 2048 MB on windows
        _prop.raw_file_size = *(__int64 *)(((char *)&(is.tellg())) + 8);
#else
        _prop.raw_file_size = is.tellg();
#endif
        is.seekg( 0, is.beg );

        std::vector<char> raw_timestep;
        raw_timestep.resize(_prop.raw_file_size);

        // read data as a block:
        is.read(raw_timestep.data(), _prop.raw_file_size);
        _raw_data.push_back(std::move(raw_timestep));

        if (!is)
            throw std::runtime_error("Error reading " + raw_file_name);
        is.close();
    }
    else
    {
        throw std::runtime_error("Could not open " + raw_file_name);
    }

    // if format was not specified in .dat file, try to calculate from
    // file size and volume resolution
    if (_prop.format.empty() && !_raw_data.empty())
    {
        unsigned int bytes = _raw_data.at(0).size() / (static_cast<size_t>(_prop.volume_res[0]) *
                                                       static_cast<size_t>(_prop.volume_res[1]) *
                                                       static_cast<size_t>(_prop.volume_res[2]));
        switch (bytes)
        {
        case 1:
            _prop.format = "UCHAR";
            std::cout << "Format determined as UCHAR." << std::endl;
            break;
        case 2:
            _prop.format = "USHORT";
            std::cout << "Format determined as USHORT." << std::endl;
            break;
        case 4:
            _prop.format = "FLOAT";
            std::cout << "Format determined as FLOAT." << std::endl;
            break;
        default: throw std::runtime_error("Could not resolve missing format specification.");
        }
    }
}
