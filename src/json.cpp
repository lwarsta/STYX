/**
 * @file json.cpp
 * @brief Implementation of JsonValue class methods and internal JSON parser.
 *
 * Implements recursive descent parsing logic for JSON structures.
 */

#include "json.h"
#include <cctype>
#include <iostream>
#include <stdexcept>

// Constructors implementation
JsonValue::JsonValue() : value(nullptr) {}
JsonValue::JsonValue(std::nullptr_t) : value(nullptr) {}
JsonValue::JsonValue(bool b) : value(b) {}
JsonValue::JsonValue(int i) : value(i) {}
JsonValue::JsonValue(double d) : value(d) {}
JsonValue::JsonValue(const std::string& s) : value(s) {}
JsonValue::JsonValue(const char* s) : value(std::string(s)) {}
JsonValue::JsonValue(const object& o) : value(o) {}
JsonValue::JsonValue(const array& a) : value(a) {}

// Type-checking methods
bool JsonValue::is_null() const { return std::holds_alternative<std::nullptr_t>(value); }
bool JsonValue::is_bool() const { return std::holds_alternative<bool>(value); }
bool JsonValue::is_int() const { return std::holds_alternative<int>(value); }
bool JsonValue::is_double() const { return std::holds_alternative<double>(value); }
bool JsonValue::is_string() const { return std::holds_alternative<std::string>(value); }
bool JsonValue::is_object() const { return std::holds_alternative<object>(value); }
bool JsonValue::is_array() const { return std::holds_alternative<array>(value); }

// Type conversion methods (accessors)
const JsonValue::object& JsonValue::as_object() const { return std::get<object>(value); }
const JsonValue::array& JsonValue::as_array() const { return std::get<array>(value); }
const std::string& JsonValue::as_string() const { return std::get<std::string>(value); }
int JsonValue::as_int() const {
    if (std::holds_alternative<int>(value))
        return std::get<int>(value);
    else if (std::holds_alternative<double>(value))
        return static_cast<int>(std::get<double>(value));
    else
        throw std::runtime_error("Value is not an int");
}
double JsonValue::as_double() const {
    if (std::holds_alternative<double>(value))
        return std::get<double>(value);
    else if (std::holds_alternative<int>(value))
        return static_cast<double>(std::get<int>(value));
    else
        throw std::runtime_error("JsonValue is not a number convertible to double");
}
bool JsonValue::as_bool() const { return std::get<bool>(value); }

/**
 * @brief Internal recursive descent parser class for JSON strings.
 */
class JsonValue::Parser {
public:
    /**
     * @brief Construct a parser with the given input string.
     * @param s Input JSON string.
     */
    Parser(const std::string& s) : str(s), pos(0) {}

    /**
     * @brief Parses the JSON value from the current position.
     * @return Parsed JsonValue.
     */
    JsonValue parse_value();

private:
    const std::string& str; ///< Input JSON string to parse.
    size_t pos;             ///< Current position in the input string.

    /**
     * @brief Returns the current character without advancing the position.
     * @return Current character or '\0' if at end of string.
     */
    char peek() const;
    /**
     * @brief Returns the current character and advances the position.
     * @return Current character.
     */
    char get();
    /**
     * @brief Skips whitespace characters.
     */
    void skip_whitespace();

    /**
     * @brief Checks if the next characters match the given literal.
     * @param literal The string to match.
     * @return True if matched, false otherwise.
     */
    bool match(const std::string& literal);
    /**
     * @brief Parses a JSON string.
     * @return The parsed string.
     */
    std::string parse_string();
    /**
     * @brief Parses a JSON number.
     * @return A JsonValue representing the number.
     */
    JsonValue parse_number();
    /**
     * @brief Parses a JSON object.
     * @return A JsonValue::object representing the parsed object.
     */
    object parse_object();
    /**
     * @brief Parses a JSON array.
     * @return A JsonValue::array representing the parsed array.
     */
    array parse_array();
};

/**
 * @brief Returns the current character without advancing the parsing position.
 * 
 * @return char The current character in the input string, or '\0' if at the end.
 */
char JsonValue::Parser::peek() const {
    return pos < str.size() ? str[pos] : '\0';
}

/**
 * @brief Returns the current character and advances the parsing position.
 * 
 * @return char The current character from the input string.
 */
char JsonValue::Parser::get() {
    return pos < str.size() ? str[pos++] : '\0';
}

/**
 * @brief Skips whitespace characters in the input string.
 */
void JsonValue::Parser::skip_whitespace() {
    while (pos < str.size() && std::isspace(static_cast<unsigned char>(str[pos])))
        pos++;
}

/**
 * @brief Attempts to match the next characters in the input string with a given literal.
 * 
 * If the literal matches, the parser advances the position by the literal's length.
 * 
 * @param literal The string to match.
 * @return true if the literal matches; false otherwise.
 */
bool JsonValue::Parser::match(const std::string& literal) {
    skip_whitespace();
    if (str.compare(pos, literal.size(), literal) == 0) {
        pos += literal.size();
        return true;
    }
    return false;
}

/**
 * @brief Parses a JSON string.
 * 
 * Expects the current position to be at the beginning of a string (i.e. a double-quote).
 * This minimal implementation does not handle escaped characters.
 * 
 * @return std::string The parsed string.
 * @throws std::runtime_error if the string does not begin with a double-quote.
 */
std::string JsonValue::Parser::parse_string() {
    if (get() != '\"')
        throw std::runtime_error("Expected '\"' at beginning of string");
    std::string result;
    while (pos < str.size()) {
        char c = get();
        if (c == '\"')
            break;
        // (Note: Escaped characters are not handled in this minimal implementation.)
        result.push_back(c);
    }
    return result;
}

/**
 * @brief Parses a JSON number.
 * 
 * Accepts digits, decimal point, sign characters, and exponent markers.
 * First attempts to parse an integer; if that fails, parses a double.
 * 
 * @return JsonValue A JsonValue representing the parsed number.
 */
JsonValue JsonValue::Parser::parse_number() {
    skip_whitespace();
    size_t start = pos;
    // Accept digits, decimal point, sign characters, and exponent markers.
    while (pos < str.size() &&
           (std::isdigit(str[pos]) || str[pos] == '.' ||
            str[pos] == '-' || str[pos] == '+' ||
            str[pos] == 'e' || str[pos] == 'E')) {
        pos++;
    }
    std::string numStr = str.substr(start, pos - start);
    // Always parse as double.
    double doubleVal = std::stod(numStr);
    return JsonValue(doubleVal);
}

/**
 * @brief Parses a JSON object.
 * 
 * Expects the current position to be at the beginning of an object ('{').
 * Recursively parses key-value pairs and returns them in a map.
 * 
 * @return JsonValue::object A map representing the parsed JSON object.
 * @throws std::runtime_error if the object syntax is invalid.
 */
JsonValue::object JsonValue::Parser::parse_object() {
    object obj;
    // Consume '{'
    if (get() != '{')
        throw std::runtime_error("Expected '{' at beginning of object");
    skip_whitespace();
    if (peek() == '}') {
        get(); // consume '}'
        return obj;
    }
    while (true) {
        skip_whitespace();
        if (peek() != '\"')
            throw std::runtime_error("Expected string key in object");
        std::string key = parse_string();
        skip_whitespace();
        if (get() != ':')
            throw std::runtime_error("Expected ':' after key in object");
        skip_whitespace();
        JsonValue value = parse_value();
        obj[key] = value;
        skip_whitespace();
        char c = get();
        if (c == '}')
            break;
        else if (c != ',')
            throw std::runtime_error("Expected ',' or '}' in object");
    }
    return obj;
}

/**
 * @brief Parses a JSON array.
 * 
 * Expects the current position to be at the beginning of an array ('[').
 * Recursively parses the array elements.
 * 
 * @return JsonValue::array A vector representing the parsed JSON array.
 * @throws std::runtime_error if the array syntax is invalid.
 */
JsonValue::array JsonValue::Parser::parse_array() {
    array arr;
    // Consume '['
    if (get() != '[')
        throw std::runtime_error("Expected '[' at beginning of array");
    skip_whitespace();
    if (peek() == ']') {
        get(); // consume ']'
        return arr;
    }
    while (true) {
        skip_whitespace();
        JsonValue value = parse_value();
        arr.push_back(value);
        skip_whitespace();
        char c = get();
        if (c == ']')
            break;
        else if (c != ',')
            throw std::runtime_error("Expected ',' or ']' in array");
    }
    return arr;
}

/**
 * @brief Parses a JSON value from the current position.
 * 
 * Determines the type of the value by inspecting the next character and
 * calls the appropriate parsing function.
 * 
 * @return JsonValue The parsed JSON value.
 * @throws std::runtime_error if the value is unrecognized.
 */
JsonValue JsonValue::Parser::parse_value() {
    skip_whitespace();
    char c = peek();
    if (c == 'n') {
        if (match("null"))
            return JsonValue(nullptr);
        else
            throw std::runtime_error("Invalid token, expected 'null'");
    } else if (c == 't') {
        if (match("true"))
            return JsonValue(true);
        else
            throw std::runtime_error("Invalid token, expected 'true'");
    } else if (c == 'f') {
        if (match("false"))
            return JsonValue(false);
        else
            throw std::runtime_error("Invalid token, expected 'false'");
    } else if (c == '\"') {
        return JsonValue(parse_string());
    } else if (c == '-' || std::isdigit(c)) {
        return parse_number();
    } else if (c == '{') {
        return JsonValue(parse_object());
    } else if (c == '[') {
        return JsonValue(parse_array());
    }
    throw std::runtime_error("Unrecognized JSON value");
}

/**
 * @brief Parses a JSON-formatted string into a JsonValue object.
 * 
 * This static function creates a Parser instance for the given string and
 * initiates the parsing process.
 * 
 * @param json_str The JSON-formatted input string.
 * @return JsonValue The parsed JSON value.
 */
JsonValue JsonValue::parse(const std::string& json_str) {
    Parser parser(json_str);
    return parser.parse_value();
}

/**
 * @brief Recursively prints the JSON value in a human-readable format.
 *
 * If the JSON value is an object or an array, it iterates over its elements and prints them with proper indentation.
 * For primitive types (string, int, double, bool, null), it prints the value.
 *
 * @param j The JsonValue to print.
 * @param indent The current indentation level (default is 0).
 */
void JsonValue::print(const JsonValue &j, int indent) {
    std::string ind(indent, ' ');
    if (j.is_object()) {
        std::cout << ind << "{\n";
        for (const auto &kv : j.as_object()) {
            std::cout << ind << "  \"" << kv.first << "\": ";
            print(kv.second, indent + 2);
        }
        std::cout << ind << "}\n";
    } else if (j.is_array()) {
        std::cout << ind << "[\n";
        for (const auto &item : j.as_array()) {
            print(item, indent + 2);
        }
        std::cout << ind << "]\n";
    } else if (j.is_string()) {
        std::cout << ind << "\"" << j.as_string() << "\"\n";
    } else if (j.is_int()) {
        std::cout << ind << j.as_int() << "\n";
    } else if (j.is_double()) {
        std::cout << ind << j.as_double() << "\n";
    } else if (j.is_bool()) {
        std::cout << ind << (j.as_bool() ? "true" : "false") << "\n";
    } else if (j.is_null()) {
        std::cout << ind << "null\n";
    }
}