// M0.10-T03: Minimal JSON parser (header-only, no external dependency)
//
// Provides a lightweight JSON value tree and streaming parser sufficient for
// reading the case_setup_json and gap_report_json string attributes written by
// the Python HDF5 writer.
#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <cstdint>
#include <utility>

namespace cfdx {
namespace io {
namespace mini_json {

class value;

using object_t = std::map<std::string, value>;
using array_t  = std::vector<value>;

enum class value_type {
    null_t,
    boolean,
    number,
    string,
    array,
    object
};

class value {
public:
    value_type type = value_type::null_t;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    array_t array;
    object_t object;

    value() = default;
    value(value_type t) : type(t) {}
    value(const std::string& s) : type(value_type::string), string(s) {}
    value(const char* s) : type(value_type::string), string(s) {}
    value(bool b) : type(value_type::boolean), boolean(b) {}
    value(double n) : type(value_type::number), number(n) {}
    value(std::int64_t n) : type(value_type::number), number(static_cast<double>(n)) {}
    value(std::size_t n) : type(value_type::number), number(static_cast<double>(n)) {}

    value& operator=(const std::string& s) {
        type = value_type::string;
        string = s;
        return *this;
    }
    value& operator=(const char* s) {
        type = value_type::string;
        string = s;
        return *this;
    }
    value& operator=(bool b) {
        type = value_type::boolean;
        boolean = b;
        return *this;
    }
    value& operator=(double n) {
        type = value_type::number;
        number = n;
        return *this;
    }
    value& operator=(const array_t& a) {
        type = value_type::array;
        array = a;
        return *this;
    }
    value& operator=(array_t&& a) {
        type = value_type::array;
        array = std::move(a);
        return *this;
    }
    value& operator=(const object_t& o) {
        type = value_type::object;
        object = o;
        return *this;
    }
    value& operator=(object_t&& o) {
        type = value_type::object;
        object = std::move(o);
        return *this;
    }

    static value parse(const std::string& text) {
        std::size_t pos = 0;
        skip_ws(text, pos);
        value v = parse_value(text, pos);
        skip_ws(text, pos);
        if (pos != text.size()) {
            throw std::runtime_error("mini_json: trailing characters at position " + std::to_string(pos));
        }
        return v;
    }

    static value parse_safe(const std::string& text) {
        try {
            return parse(text);
        } catch (...) {
            return value(value_type::null_t);
        }
    }

    bool has(const std::string& key) const {
        return type == value_type::object && object.find(key) != object.end();
    }

    const value& at(const std::string& key) const {
        if (type != value_type::object) throw std::runtime_error("mini_json: not an object");
        auto it = object.find(key);
        if (it == object.end()) throw std::runtime_error("mini_json: key not found: " + key);
        return it->second;
    }

    const value& at(std::size_t i) const {
        if (type != value_type::array) throw std::runtime_error("mini_json: not an array");
        if (i >= array.size()) throw std::runtime_error("mini_json: array index out of range");
        return array[i];
    }

    bool is_object() const noexcept { return type == value_type::object; }
    bool is_array() const noexcept { return type == value_type::array; }
    bool is_string() const noexcept { return type == value_type::string; }
    bool is_number() const noexcept { return type == value_type::number; }
    bool is_boolean() const noexcept { return type == value_type::boolean; }
    bool is_null() const noexcept { return type == value_type::null_t; }

    std::string as_string() const {
        if (type != value_type::string) throw std::runtime_error("mini_json: not a string");
        return string;
    }

    double as_number() const {
        if (type != value_type::number) throw std::runtime_error("mini_json: not a number");
        return number;
    }

    std::int64_t as_int() const {
        return static_cast<std::int64_t>(as_number());
    }

    bool as_boolean() const {
        if (type != value_type::boolean) throw std::runtime_error("mini_json: not a boolean");
        return boolean;
    }

    std::string serialize() const {
        std::ostringstream ss;
        write_to(ss);
        return ss.str();
    }

    void write_to(std::ostringstream& ss) const {
        switch (type) {
            case value_type::null_t:   ss << "null"; break;
            case value_type::boolean:  ss << (boolean ? "true" : "false"); break;
            case value_type::number:   {
                if (number == static_cast<double>(static_cast<std::int64_t>(number))) {
                    ss << static_cast<std::int64_t>(number);
                } else {
                    ss << number;
                }
                break;
            }
            case value_type::string:   ss << '"' << escape_string(string) << '"'; break;
            case value_type::array: {
                ss << '[';
                for (std::size_t i = 0; i < array.size(); ++i) {
                    if (i > 0) ss << ',';
                    array[i].write_to(ss);
                }
                ss << ']';
                break;
            }
            case value_type::object: {
                ss << '{';
                bool first = true;
                for (const auto& [k, v] : object) {
                    if (!first) ss << ',';
                    first = false;
                    ss << '"' << escape_string(k) << "\":";
                    v.write_to(ss);
                }
                ss << '}';
                break;
            }
        }
    }

    static std::string escape_string(const std::string& s) {
        std::string out;
        out.reserve(s.size() + 8);
        for (char c : s) {
            switch (c) {
                case '"':  out += "\\\""; break;
                case '\\': out += "\\\\"; break;
                case '\n': out += "\\n";  break;
                case '\r': out += "\\r";  break;
                case '\t': out += "\\t";  break;
                default:   out += c;      break;
            }
        }
        return out;
    }

private:
    static void skip_ws(const std::string& text, std::size_t& pos) {
        while (pos < text.size()) {
            char c = text[pos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                ++pos;
            } else {
                break;
            }
        }
    }

    static value parse_value(const std::string& text, std::size_t& pos) {
        skip_ws(text, pos);
        if (pos >= text.size()) throw std::runtime_error("mini_json: unexpected end of input");

        char c = text[pos];
        if (c == '{') return parse_object(text, pos);
        if (c == '[') return parse_array(text, pos);
        if (c == '"') return parse_string(text, pos);
        if (c == 't' || c == 'f') return parse_boolean(text, pos);
        if (c == 'n') return parse_null(text, pos);
        if (c == '-' || (c >= '0' && c <= '9')) return parse_number(text, pos);
        throw std::runtime_error("mini_json: unexpected character '" + std::string(1, c) + "' at position " + std::to_string(pos));
    }

    static value parse_object(const std::string& text, std::size_t& pos) {
        value v(value_type::object);
        ++pos; // skip '{'
        skip_ws(text, pos);
        if (pos < text.size() && text[pos] == '}') {
            ++pos;
            return v;
        }
        while (true) {
            skip_ws(text, pos);
            std::string key = parse_string_raw(text, pos);
            skip_ws(text, pos);
            if (pos >= text.size() || text[pos] != ':') {
                throw std::runtime_error("mini_json: expected ':' at position " + std::to_string(pos));
            }
            ++pos; // skip ':'
            value val = parse_value(text, pos);
            v.object[key] = std::move(val);
            skip_ws(text, pos);
            if (pos < text.size() && text[pos] == ',') {
                ++pos;
                continue;
            }
            if (pos < text.size() && text[pos] == '}') {
                ++pos;
                break;
            }
            throw std::runtime_error("mini_json: expected ',' or '}' at position " + std::to_string(pos));
        }
        return v;
    }

    static value parse_array(const std::string& text, std::size_t& pos) {
        value v(value_type::array);
        ++pos; // skip '['
        skip_ws(text, pos);
        if (pos < text.size() && text[pos] == ']') {
            ++pos;
            return v;
        }
        while (true) {
            value elem = parse_value(text, pos);
            v.array.push_back(std::move(elem));
            skip_ws(text, pos);
            if (pos < text.size() && text[pos] == ',') {
                ++pos;
                continue;
            }
            if (pos < text.size() && text[pos] == ']') {
                ++pos;
                break;
            }
            throw std::runtime_error("mini_json: expected ',' or ']' at position " + std::to_string(pos));
        }
        return v;
    }

    static std::string parse_string_raw(const std::string& text, std::size_t& pos) {
        skip_ws(text, pos);
        if (pos >= text.size() || text[pos] != '"') {
            throw std::runtime_error("mini_json: expected string at position " + std::to_string(pos));
        }
        ++pos; // skip opening quote
        std::string result;
        while (pos < text.size()) {
            char c = text[pos];
            if (c == '"') {
                ++pos;
                return result;
            }
            if (c == '\\') {
                ++pos;
                if (pos >= text.size()) throw std::runtime_error("mini_json: truncated escape sequence");
                char esc = text[pos];
                switch (esc) {
                    case '"':  result += '"';  break;
                    case '\\': result += '\\'; break;
                    case '/':  result += '/';  break;
                    case 'n':  result += '\n'; break;
                    case 'r':  result += '\r'; break;
                    case 't':  result += '\t'; break;
                    case 'b':  result += '\b'; break;
                    case 'f':  result += '\f'; break;
                    case 'u': {
                        if (pos + 4 >= text.size()) throw std::runtime_error("mini_json: invalid unicode escape");
                        std::string hex = text.substr(pos + 1, 4);
                        unsigned int codepoint = static_cast<unsigned int>(std::stoul(hex, nullptr, 16));
                        if (codepoint < 0x80) {
                            result += static_cast<char>(codepoint);
                        } else if (codepoint < 0x800) {
                            result += static_cast<char>(0xC0 | (codepoint >> 6));
                            result += static_cast<char>(0x80 | (codepoint & 0x3F));
                        } else {
                            result += static_cast<char>(0xE0 | (codepoint >> 12));
                            result += static_cast<char>(0x80 | ((codepoint >> 6) & 0x3F));
                            result += static_cast<char>(0x80 | (codepoint & 0x3F));
                        }
                        pos += 4;
                        break;
                    }
                    default:
                        result += esc;
                        break;
                }
                ++pos;
            } else {
                result += c;
                ++pos;
            }
        }
        throw std::runtime_error("mini_json: unterminated string");
    }

    static value parse_string(const std::string& text, std::size_t& pos) {
        value v(value_type::string);
        v.string = parse_string_raw(text, pos);
        return v;
    }

    static value parse_boolean(const std::string& text, std::size_t& pos) {
        value v(value_type::boolean);
        if (text.compare(pos, 4, "true") == 0) {
            v.boolean = true;
            pos += 4;
        } else if (text.compare(pos, 5, "false") == 0) {
            v.boolean = false;
            pos += 5;
        } else {
            throw std::runtime_error("mini_json: invalid boolean at position " + std::to_string(pos));
        }
        return v;
    }

    static value parse_null(const std::string& text, std::size_t& pos) {
        if (text.compare(pos, 4, "null") == 0) {
            pos += 4;
            return value(value_type::null_t);
        }
        throw std::runtime_error("mini_json: invalid null at position " + std::to_string(pos));
    }

    static value parse_number(const std::string& text, std::size_t& pos) {
        std::size_t start = pos;
        if (pos < text.size() && text[pos] == '-') ++pos;
        while (pos < text.size() && text[pos] >= '0' && text[pos] <= '9') ++pos;
        if (pos < text.size() && text[pos] == '.') {
            ++pos;
            while (pos < text.size() && text[pos] >= '0' && text[pos] <= '9') ++pos;
        }
        if (pos < text.size() && (text[pos] == 'e' || text[pos] == 'E')) {
            ++pos;
            if (pos < text.size() && (text[pos] == '+' || text[pos] == '-')) ++pos;
            while (pos < text.size() && text[pos] >= '0' && text[pos] <= '9') ++pos;
        }
        value v(value_type::number);
        try {
            v.number = std::stod(text.substr(start, pos - start));
        } catch (...) {
            v.number = 0.0;
        }
        return v;
    }
};

// Convenience: safely access nested value, returning nullptr if not found
inline const value* find(const value& v, const std::string& key) {
    if (v.is_object() && v.has(key)) {
        return &v.at(key);
    }
    return nullptr;
}

}  // namespace mini_json
}  // namespace io
}  // namespace cfdx