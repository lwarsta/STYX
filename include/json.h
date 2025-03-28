/**
 * @file json.h
 * @brief Single-class JSON parser and interpreter supporting nested objects, arrays, strings, integers, doubles, and booleans.
 *
 * @author
 * @date 2024
 */

#ifndef JSON_H
#define JSON_H

#include <variant>
#include <string>
#include <map>
#include <vector>

/**
 * @brief Class representing a JSON value.
 *
 * Supports JSON null, boolean, integer, double, string, object, and array types.
 */
class JsonValue {
public:
    /// JSON object representation (key-value pairs)
    using object = std::map<std::string, JsonValue>;

    /// JSON array representation (ordered values)
    using array = std::vector<JsonValue>;

    /// Variant type holding the actual JSON value
    std::variant<std::nullptr_t, bool, int, double, std::string, object, array> value;

    /// Constructors for various JSON value types
    JsonValue();
    JsonValue(std::nullptr_t);
    JsonValue(bool b);
    JsonValue(int i);
    JsonValue(double d);
    JsonValue(const std::string& s);
    JsonValue(const char* s);
    JsonValue(const object& o);
    JsonValue(const array& a);

    /// Type-checking functions
    bool is_null() const;
    bool is_bool() const;
    bool is_int() const;
    bool is_double() const;
    bool is_string() const;
    bool is_object() const;
    bool is_array() const;

    /// Type conversion functions (accessors)
    const object& as_object() const;
    const array& as_array() const;
    const std::string& as_string() const;
    int as_int() const;
    double as_double() const;
    bool as_bool() const;

    /**
     * @brief Parse JSON-formatted string into a JsonValue object.
     *
     * @param json_str JSON-formatted input string.
     * @return Parsed JsonValue object.
     */
    static JsonValue parse(const std::string& json_str);
	
	/**
	 * @brief Recursively prints the JSON value in a human-readable format.
	 *
	 * If the JSON value is an object or an array, it iterates over its elements and prints them with proper indentation.
	 * For primitive types (string, int, double, bool, null), it prints the value.
	 *
	 * @param j The JsonValue to print.
	 * @param indent The current indentation level (default is 0).
	 */
	void print(const JsonValue &j, int indent = 0);

private:
    /// Forward declaration of the parser class
    class Parser;
};

#endif