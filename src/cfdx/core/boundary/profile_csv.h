#pragma once

#include "cfdx/core/boundary/value_provider.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

inline std::shared_ptr<const ValueProvider> read_spatial_profile_csv(
    const std::string& path, int axis = 1)
{
    std::ifstream input(path);
    if (!input.is_open())
        throw std::invalid_argument("inlet profile: cannot open CSV '" + path + "'");

    const auto trim = [](std::string value) {
        const auto first = std::find_if_not(value.begin(), value.end(),
            [](unsigned char ch) { return std::isspace(ch) != 0; });
        const auto last = std::find_if_not(value.rbegin(), value.rend(),
            [](unsigned char ch) { return std::isspace(ch) != 0; }).base();
        if (first >= last) return std::string{};
        return std::string(first, last);
    };
    const auto normalized = [&](std::string value) {
        value = trim(std::move(value));
        std::transform(value.begin(), value.end(), value.begin(),
            [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        return value;
    };
    const auto parse_number = [&](const std::string& token,
                                  std::size_t line_number) {
        const std::string value = trim(token);
        std::size_t consumed = 0;
        double parsed = 0.0;
        try {
            parsed = std::stod(value, &consumed);
        } catch (const std::exception&) {
            throw std::invalid_argument(
                "inlet profile: invalid numeric field at line " +
                std::to_string(line_number));
        }
        while (consumed < value.size() &&
               std::isspace(static_cast<unsigned char>(value[consumed])) != 0)
            ++consumed;
        if (consumed != value.size() || !std::isfinite(parsed))
            throw std::invalid_argument(
                "inlet profile: malformed or non-finite numeric field at line " +
                std::to_string(line_number));
        return parsed;
    };

    std::vector<double> coordinates, values;
    std::string line;
    std::size_t line_number = 0;
    bool first_record = true;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string record = trim(line);
        if (record.empty() || record.front() == '#') continue;

        std::istringstream fields(record);
        std::string coordinate_token, value_token, extra;
        if (!std::getline(fields, coordinate_token, ',') ||
            !std::getline(fields, value_token, ',') ||
            std::getline(fields, extra, ','))
            throw std::invalid_argument(
                "inlet profile: expected exactly two CSV columns at line " +
                std::to_string(line_number));

        if (first_record) {
            const auto header_x = normalized(coordinate_token);
            const auto header_y = normalized(value_token);
            const bool recognized_header =
                (header_x == "coord" || header_x == "coordinate" ||
                 header_x == "x" || header_x == "y" || header_x == "z") &&
                (header_y == "value" || header_y == "u" ||
                 header_y == "velocity");
            if (recognized_header) {
                first_record = false;
                continue;
            }
        }
        first_record = false;
        coordinates.push_back(parse_number(coordinate_token, line_number));
        values.push_back(parse_number(value_token, line_number));
    }

    if (coordinates.size() < 2)
        throw std::invalid_argument(
            "inlet profile: CSV '" + path + "' needs at least two numeric rows");
    return std::make_shared<const TableOverSpaceValueProvider>(
        std::move(coordinates), std::move(values), axis);
}

} // namespace cfdx::core
