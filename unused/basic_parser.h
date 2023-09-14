#ifndef basic_parse_hh
#define basic_parse_hh


#include <string>
#include <sstream>
#include <map>
#include <array>
#include <vector>
#include <fstream>
#include <cctype>
#include <iostream>
#include <filesystem>
#include <algorithm>


/**
 * A barebones parameter parser for scientific simulations.
 * 
 * This is intended as a drop-in replacement for long strings of command line
 * arguments that lose their meaning and become confusing.
 * 
 * The main tools are `declare` and `save_file`:
 * 
 * Name your temperature with
 * 
 *      basic_parser p;
 * 
 *      double T; 
 *      p.declare("Temperature (K)", T);
 * 
 *      then 
 *      p.from_file('input.txt')
 * 
 * where 'input.txt' reads
 * 
 *      # this is a comment
 *      Temperature (K): 420
 * 
 * 
 * This basically implents a subset of the TOML standard, but with only top-level directories.
 * 
 * Only top level structure is supported, for more complicated inputs
 * you should just use a library for TOML/YAML/JSON files.
 * 
 */

/**
 * @brief Returns current date and time formatted as YYYY-MM-DD.HH:mm:ss
 * @return current date/time, format is YYYY-MM-DD.HH:mm:ss
 */
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%dT%X", &tstruct);

    return buf;
}

/**
 * @brief Gets the current git commit (prints a warning if unable)
 * 
 * @return const std::string current git commit hash
 */
const std::string currentCommit(){
    const char* git_log = "git --no-pager log --pretty=format:'%h' -n 1";

    std::array<char, 128> buffer;
    std::string result = "#";
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(git_log, "r"), pclose);
    if (!pipe) {
        result += "xxxxxxxx";
        std::cerr<<"[ WARN ] Could not determine git hash, padding with 'xxxxxxxx' \n";
    } else {
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
            result += buffer.data();
        }
    }
    return result;    
}

/**
 * @brief Like a python type, just a label
 * 
 */
enum class paramtype {
    None, Int, UInt, Float, String, Bool, Enum
};

std::string to_string(paramtype p){
    switch (p)
    {
    case paramtype::Int:
        return "Integer";
    case paramtype::UInt:
        return "Unsigned Integer";
    case paramtype::Float:
        return "Float";
    case paramtype::String:
        return "String";
    case paramtype::Bool:
        return "Boolean";
    case paramtype::Enum:
        return "Enum";
    default:
        return "NoneType";
    }
}


/**
 * @brief Basic parameter parser
 * This works in the procedural paradigm - it attaches references from outside to string handles,
 * which must be named at compile time. The only state stored by this class is the set of flags 
 * indicating whether or not basic_parser have been initialised.
 * 
 * @tparam bp_int_t type to use for signed int
 * @tparam bp_uint_t type to use for unsigned int
 * @tparam bp_float_t type to use for floating points
 * 
 */
template <typename bp_int_t=int64_t, typename bp_uint_t=uint64_t, typename bp_float_t=double>
class basic_parser {
    template <typename it, typename uit, typename ft>
    friend std::ostream& operator<<(std::ostream& os, basic_parser<it, uit, ft> p);
    template <typename it, typename uit, typename ft>
    friend std::istream& operator>>(std::istream& is, basic_parser<it, uit, ft> p);
public:
    /**
     * @brief Construct a new basic parser object
     * 
     * @param header Line to open output files with under save_input.
     * @param comment_char Lines starting with this character in infiles are ignored.
     * @param output_char Token delimiting parameter assignment in filenames
     */
    basic_parser(const std::string& header="", char comment_char='#', char output_char='%') : 
        header(header), comment_char(comment_char), output_delimiter(output_char) {
        };

    /**
     * @brief Outputs current parser state into human readable text
     * 
     * @param os std::ostream to dump data into
     * @param delimiter What should go between the data handle and the data, e.g. `handle::data`
     */
    void into_stream(std::ostream& os, const char* delimiter="=");

    /**
     * @brief Parses data from stream
     * 
     * @param is std::istream to read data from
     * @param delimiter What should go between the data handle and the data, e.g. `handle::data`
     */
    void from_stream(std::istream& is, const char* delimiter="=");

    /**
     * @brief Reads an input file, setting all of the known variables using it
     * 
     * @param fname Filename (C string)
     * @param delimiter Delimiter (C string)
     */
    void from_file(const std::filesystem::path& fname, const char* delimiter="=");
    void from_file(const char* fname, const char* delimiter="="){
        from_file(std::filesystem::path(fname), delimiter);
    }
    void from_file(const std::string& fname, const char* delimiter="="){
        from_file(std::filesystem::path(fname), delimiter);
    }

    /**
     * @brief Interprets command line options in the form --handle=55
     * 
     * @param argc The end of the argv array
     * @param argv it's argv
     * @param start=1 Expects argv[start] to be "--handle1".
     * 
     * @return a pretty-printed string of which params were set
     */
    std::string from_argv(int argc, const char** argv, int start=1);

    /**
     * @brief Reads an input file, setting all of the known variables using it
     * 
     * @param fname Filename (C string)
     * @param delimiter Delimiter (C string)
     */
    void into_file(const char* fname, const char* delimiter="=");

    void into_file(const std::string& fname, const char* delimiter="="){
        into_file(fname.c_str(), delimiter);
    }


    /**
     * @brief Throws an exception if parameters have not been set
     * 
     */
    void assert_initialised();

    // /**
    //  * @brief Only passes silently if one of the following conditions are met:
    //  *   1. All parameters in v1 are initialised AND no parameters of v2 are initialised.
    //  *   2. All parameters in v2 are initialised AND no parameters of v1 are initialised.
    //  *   All parameters MUST have been declared to make this valid.
    //  * @param v1 
    //  * @param v2 
    //  */
    // void assert_exclusive(const std::vector<std::string>& v1, const std::vector<std::string>& v2);

    /**
     * @brief Associates the C string 'handle' with the variable x.
     * 
     * @param handle The C string to search for in the infile.
     * @param x The varaible that the parser should store the result in (passed by reference)
     */
    void declare(const std::string& handle, bp_int_t* x){
        assert_unique(handle);
        ints[handle] = x;
        index[handle] = paramtype::Int;
        initialised[handle]=false;
    }
    /// @overload
    void declare(const std::string& handle, bp_uint_t* x){
        assert_unique(handle);
        uints[handle] = x;
        index[handle] = paramtype::UInt;
        initialised[handle]=false;
    }
    /// @overload
    void declare(const std::string& handle, bp_float_t* x){
        assert_unique(handle);
        floats[handle] = x;
        index[handle] = paramtype::Float;
        initialised[handle]=false;
    }
    /// @overload
    void declare(const std::string& handle, bool* b){
        assert_unique(handle);
        bools[handle] = b;
        index[handle] = paramtype::Bool;
        initialised[handle]=false;
    }
    /// @overload
    void declare(const std::string& handle, std::string* s){
        assert_unique(handle);
        strings[handle] = s;
        index[handle] = paramtype::String;
        initialised[handle]=false;
    }
    // TODO file/dir types with validation

    /**
     * @brief Declares an enumerated type
     * 
     * @tparam I  Integer castable enum type.
     * @param handle The (human-readable) name grepped for in the input file.
     * @param e A pointer to the enum.
     * @param lexicon An initialiser list, of the form {{"str1", enum_type::value1}, {"str2", enum_type::value2}, ...}
     *  This list is used to parse values of enums
     */
    template<typename I>
    void declare_enum(const std::string& handle, I* e,  std::initializer_list<std::pair<const std::string, I>> lexicon ){
        assert_unique(handle);
        index[handle] = paramtype::Enum;
        initialised[handle]=false;
        enums[handle] = reinterpret_cast<int*>(e);
        // enum_maps[handle] = std::map<const std::string, int>();
        for (auto& [opt, value] : lexicon){
            enum_maps[handle][opt] = static_cast<int>(value);
        }
    }

    template<typename I>
    void declare_optional_enum(const std::string& handle, I* e,  std::initializer_list<std::pair<const std::string, I>> lexicon, I default_e ){
        declare_enum(handle, e, lexicon);
        *e = default_e;
        initialised[handle]=true;
        
    }

    std::string enum_value(const std::string& handle, int idx){
        for (const auto& [key, value] : enum_maps[handle]) {
            if (idx == value) return key;
        }
        return "????";
    }
    

    /**
     * @brief Like declare(handle, x), with a default parameter.
     * 
     * @param handle The C string to search for in the infile.
     * @param x The variable that the parser should store the result in (passed by reference).
     * @param default_x The default value of x.
     * 
     * @see void declare(const std::string& handle, bp_int_t& x)
     */
    void declare_optional(const std::string& handle, bp_int_t* x, bp_int_t default_x){
        declare(handle, x);
        *x = default_x;
        initialised[handle]=true;
    }
    /// @overload
    void declare_optional(const std::string& handle, bp_uint_t* x, bp_uint_t default_x){
        declare(handle, x);
        *x = default_x;
        initialised[handle]=true;
    }
    /// @overload
    void declare_optional(const std::string& handle, bp_float_t* x, bp_float_t default_x){
        declare(handle, x);
        *x = default_x;
        initialised[handle]=true;
    }
    /// @overload
    void declare_optional(const std::string& handle, bool* b, bool default_b=false){
        declare(handle, b);
        *b = default_b;
        initialised[handle]=true;
    }
    /// @overload
    void declare_optional(const std::string& handle, std::string* s, const std::string& default_s){
        declare(handle, s);
        *s = default_s;
        initialised[handle]=true;
    }
    


    /**
     * @brief Lists all parameter handles.
     * 
     * @return A nicely formatted string of all names expected from the input file.
     */
    std::string cparam_names(const char* delimiter="="){
        std::stringstream s;
        for (auto& [handle, v] : ints) {
            s <<  handle << delimiter << "\t[ int ] ";
            if (initialised[handle]) s<<*v;
            s << "\n" ;
        }
        for (auto& [handle, v] : uints) {
            s <<  handle << delimiter << "\t[uint ] ";
            if (initialised[handle]) s<<*v;
            s << "\n" ;
        }
        for (auto& [handle, v] : floats) {
            s <<  handle << delimiter << "\t[float] ";
            if (initialised[handle]) s<<*v;
            s << "\n" ;
        }
        for (auto& [handle, v] : bools) {
            s <<  handle << delimiter << "\t[bool ] ";
            if (initialised[handle]) s<<*v;
            s << "\n" ;
        }
        for (auto& [handle, v] : strings) {
            s <<  handle << delimiter << "\t[ str ] ";
            if (initialised[handle]) s<<*v;
            s << "\n" ;
        }
        return s.str();
    }

    /**
     * @brief Get the outfile prefix 
     * 
     * @return std::string formed by concatenating the extensionless parameter file name with %-terminated overrides
     */
    std::string get_outfile_prefix(bool run_date=true, bool gh_hash=true){
        std::string ofp;
        if (run_date) {ofp += currentDateTime(); }
        if (gh_hash) { ofp += currentCommit(); }
        // todo find a good way to include build date  { ofp += __DATE__; }
        
        for (int j=0; j<this->loaded_files.size(); j++){
            ofp += this->loaded_files[j];
            if (j != 0){
                ofp += "+";
            } 
        }
        for (auto& [p, v] : this->argv_overrides){
            ofp += output_delimiter + p + "=" + v;
        }
        return ofp;
    }

private:
    // the pointers
    std::map<const std::string,  bp_int_t*     > ints;
    std::map<const std::string,  bp_uint_t*    > uints;
    std::map<const std::string,  bp_float_t*   > floats;
    std::map<const std::string,  std::string*  > strings;
    std::map<const std::string,  bool*         > bools;
    std::map<const std::string,  int*          > enums;
    // std::map<const std::string,  std::filesystem::path* > files; 
    // std::map<const std::string,  std::filesystem::path* > folders;

    // the conversions string to int for enum purposes
    std::map<const std::string, std::map<std::string, int> > enum_maps;


    

    const std::string header;
    std::vector<std::string> loaded_files;
    std::vector<std::pair<std::string, std::string> > argv_overrides;
    const char comment_char;
    const char output_delimiter;

    bool set_value(const std::string& handle, const std::string& value);

    void assert_unique(const std::string& handle){
        // Make sure we aren't doubling up
        if (index.find(handle)!= index.end()){
            std::cerr<<"Duplicate Declaration found!"<<std::endl;
            throw std::range_error("Duplicate entry");
        }
    }

    std::map<const std::string, bool> initialised; /// Keeps track of whether the variable has been set or not

    std::map<const std::string, paramtype> index; /// Stores the order in which handles were passed for consistency
};

std::string strip(const std::string& str)
{
    std::string s(str);
    s.erase(0,s.find_first_not_of(" \t\n\r\f\v"));
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
    return s;
}

template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
void basic_parser<bp_int_t, bp_uint_t, bp_float_t>::from_stream(std::istream& is, const char* delimiter){
    size_t lineno = 0;
    for (std::string line; std::getline(is, line); ){
        lineno++;
        // remove all whitespace from the line
        line = strip(line);
        // Ignore comments and blank lines
        if ( line[0] == comment_char || line.size() == 0) continue;
        
        // Look for delimiters
        size_t idx = line.find(delimiter);

        // If we could not find it, raise an exception
        if (idx == std::string::npos){
            fprintf(stderr, "Invalid line found:\n%03lu | %s", lineno, line.c_str()); fflush(stderr);
            throw std::runtime_error("Invalid file format");
        }

        std::string start = strip(line.substr(0,idx));
        std::string end = strip(line.substr(idx+1));

        try
        {
            if (!set_value(start, end)){
                fprintf(stderr, "Invalid specification on line %03lu:\n '%s'\n", lineno, line.c_str()); 
                fprintf(stderr, "Expected basic_parser format:\n%s\n",cparam_names().c_str());
                fflush(stderr);
                throw std::runtime_error("Invalid specification");
            }
        }
        catch(const std::exception& e)
        {   
            std::cerr<< "Bad input file: " << e.what() << "\n"; 
            std::cerr<< "Expected basic_parser format:\n"<< cparam_names() << std::endl;
            throw e;
            
        }
       
    }
}

bool stobool(const std::string& s){
    std::string ss(s);
    for (char& c : ss){
        c = std::tolower(c);
    }
    if (s == "true"){
        return true;
    } else if (s=="false") {
        return false;
    } else {
        throw std::runtime_error("No valid conversion from string to bool");
    }
}

template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
bool basic_parser<bp_int_t, bp_uint_t, bp_float_t>::set_value(const std::string& handle, const std::string& value){
     // see if we recognise the index
        switch (index[handle])
        {
        case paramtype::Int:
            *ints[handle] = stoll(value);
            break;
        case paramtype::UInt:
            *uints[handle] = stoull(value);
            break;
        case paramtype::Float:
            if ( (value.front() == '"' && value.back() == '"') 
                || (value.front() == '\'' && value.back() == '\'') )
            {
                    *floats[handle] = stod(value.substr(1,value.length()-2));
            } else {
                *floats[handle] = stod(value);
            }
            break;
        case paramtype::Bool:
            *bools[handle] = stobool(value);
            break;
        case paramtype::String:
            if ( (value.front() == '"' && value.back() == '"') 
                || (value.front() == '\'' && value.back() == '\'') )
            {
                *strings[handle] = std::string_view(value).substr(1,value.length()-2);
            } else {
                *strings[handle] = strip(value);
            }
            break;
        case paramtype::Enum:
            try {
                *enums[handle] = enum_maps[handle].at(value);
            } catch (std::out_of_range &e ){
                throw std::out_of_range("Unexpected value for enum " + handle + ": " + value );
            }
            break;
        default:
            throw std::runtime_error("Unknown type");
            return false;
        }
        initialised[handle] = true;
        return true;
}


template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
void basic_parser<bp_int_t, bp_uint_t, bp_float_t>::into_stream(std::ostream& os, const char* delimiter){
    std::stringstream ss(header);
    os << comment_char << " --- OUTPUT FILE --- \n";
    for (std::string line; std::getline(ss, line); ){
        os<<comment_char<<line<<"\n";
    }

    for (auto& [handle, type] : index){
        os << handle <<delimiter;
        switch (type)
        {
        case paramtype::Int:
            os << *ints[handle]<<"\n";
            break;
        case paramtype::UInt:
            os << *uints[handle]<<"\n";
            break;
        case paramtype::Float:
            os << *floats[handle]<<"\n";
            break;
        case paramtype::String:
            os <<'"'<< *strings[handle]<<"\"\n";
            break;
        case paramtype::Bool:
            os << (*bools[handle]? "\"true\"" : "\"false\"")<<"\n";
            break;
        case paramtype::Enum:
            os << this->enum_value(handle, *enums[handle])<<"\n";
        default:
            os << handle <<delimiter<< "unknown! \n";
            break;
        }
    }
}

template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
std::ostream& operator<<(std::ostream& os, basic_parser<bp_int_t, bp_uint_t, bp_float_t> p){
    p.into_stream(os);
    return os;
}

template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
std::istream& operator>>(std::istream& is, basic_parser<bp_int_t, bp_uint_t, bp_float_t> p){
    p.from_stream(is);
    return is;
}

template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
void basic_parser<bp_int_t, bp_uint_t, bp_float_t>::into_file(const char* fname, const char* delimiter){
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        fprintf(stderr, "Error opening file %s\n", fname);
        throw std::runtime_error("Cannot open file");
    }
    try{
        into_stream(ofs, delimiter);
        ofs.close();
    } catch  (const char* e) {
        fprintf(stderr, "Error writing file: %s\n",e);
        throw std::runtime_error("Cannot write to file");
        ofs.close();
    }
}
    

template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
void basic_parser<bp_int_t, bp_uint_t, bp_float_t>::from_file(const std::filesystem::path& fname, const char* delimiter){
    std::ifstream ifs(fname);
    if (!ifs.is_open()) {
        fprintf(stderr, "Error opening file!\n");
    }
    try {
        from_stream(ifs, delimiter);
        ifs.close();
    } catch  (const char* e) {
        fprintf(stderr, "Error reading file: %s\n",e);
        ifs.close();
        throw e;
    }
    // Remember this name, boy
    this->loaded_files.push_back(fname.filename().replace_extension());
}


/**
 * @brief Parses command line arguments, starting from argv[start]
 * 
 * @tparam bp_int_t Integer Type
 * @tparam bp_uint_t Unsigned Integer Type
 * @tparam bp_float_t Float Type
 * @param argc Number of nonzero command line arguments (including the ones skipped by start)
 * @param argv Argument vector
 * @param start First index to check
 * @return std::string Command line overrides in format %arg1=val1%arg2=val2 etc.
 */
template<typename bp_int_t, typename bp_uint_t, typename bp_float_t>
std::string basic_parser<bp_int_t, bp_uint_t, bp_float_t>::from_argv(int argc, const char** argv, int start)
{
    std::stringstream retval;

    for (int i=start; i<argc; i++)
    {
        std::string s(argv[i]);
        if (s[0] != '-' || s[1] != '-')
        {
            throw std::runtime_error("Bad keyword argument: handles must be prefixed with '--'");
        }
        s.erase(0,2);

        // Look for delimiters
        size_t idx = s.find('=');
        // If we could not find it, raise an exception
        if (idx == std::string::npos){
            std::cerr << "Invalid kwarg: "<< s <<std::endl;
            throw std::runtime_error("Invalid kwarg");
        }
        
        std::string start = strip(s.substr(0,idx));
        std::string end = strip(s.substr(idx+1));

        try
        {
            bool success = set_value(start, end);
            if (success)
            {
                retval << output_delimiter << start << "=" << end;
                this->argv_overrides.push_back(std::pair(start, end));
            }
            else
            {
                std::cerr<< "Invalid kwarg at position" << i<<": "<<s<<"\n";
                throw std::runtime_error("Invalid kwarg");
            }
        }
        catch(const std::exception& e)
        {   
            std::cerr<< "Bad kwarg at position "<<i<<": " << e.what() << "\n"; 
            throw std::runtime_error("Invalid kwarg");
        }
    }
    return retval.str();
}


template <typename bp_int_t, typename bp_uint_t, typename bp_float_t>
void basic_parser<bp_int_t, bp_uint_t, bp_float_t>::assert_initialised() {
    bool giveup=false;
    for (const auto& [k, x] : initialised){
        if (x == false){
            fprintf(stderr, "Uninitialied %s: %s\n",to_string(index[k]).c_str(),k.c_str());
            giveup=true;
        }
    }
    if (giveup){
        fprintf(stderr, "Expected format:\n%s",cparam_names().c_str());
        throw std::runtime_error("Mandatory variables missing from infile");
    }
}




#endif