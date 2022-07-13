#ifndef __READINPUT_H_CMC__
#define __READINPUT_H_CMC__
#include <algorithm>
#include <fstream>
#include <vector>
#include <cctype>
#include <sstream>
#include <cstdlib>
#include <cxxabi.h>
#include "StringUtility.h"
using namespace std;

ifstream open_file (const string& file)
{
    ifstream ifs (file);
    if (!ifs)
    {
        cout << "Cannot open " << file << endl;
        throw;
    }
    return ifs;
}

int count_lines (const string& file)
{
    int number_of_lines = 0;
    auto ifs = open_file (file);
    std::string line;
    while (std::getline (ifs, line))
        ++number_of_lines;
    ifs.close();
    return number_of_lines;
}

template <typename T>
vector<vector<T>> readtxt (const string& file, int skipline=0)
{
    int N = count_lines (file);
    vector<vector<T>> re;
    re.reserve (N);

    auto ifs = open_file (file);
    string line;
    for(int i = 0; i < skipline; i++)
        getline(ifs, line);
    while (getline(ifs, line))
    {
        re.push_back (split_str<T> (line));
    }
    ifs.close();
    return re;
}

template<typename T>
std::string type_name()
{
    int status;
    std::string tname = typeid(T).name();
    char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if(status == 0) {
        tname = demangled_name;
        std::free(demangled_name);
    }   
    return tname;
}

template <typename T>
vector<T> read_vector (ifstream& ifs, const string& key, int skip=2)
{
    // Go to the key-word line
    string line;
    while (getline(ifs, line))
    {
        vector<string> words = split_str<string> (line);
        if (std::find (words.begin(), words.end(), key) != words.end())
            break;
    }

    // Read data
    vector<T> re = split_str<T> (line, skip);
    //if (re.size() == 0) {
    //    cout << "ReadInput.h: read_vector: " << key << " not found" << endl;
    //    cout << "             " << line << endl;
    //    throw;
    //}
    return re;
}

template <typename T>
vector<T> read_vector (const string& fname, const string& key, int skip=2)
{
    ifstream ifs (fname);
    return read_vector<T> (ifs, key, skip);
}

bool fstream_goto (ifstream& ifs, const string& key, int skipline=0, string end="")
{
    // Go to the key-word line
    string line;
    while (getline(ifs, line))
    {
        vector<string> words = split_str<string> (line);

        if (std::find (words.begin(), words.end(), key) != words.end())
            break;
        if (end != "")
            if (std::find (words.begin(), words.end(), end) != words.end())
                return false;
        //auto n = line.find (key);
        //if (n != string::npos) break;
    }

    if (!ifs) return false;

    for(int i = 0; i < skipline; i++)
        getline(ifs, line);
    return true;
}

bool has_keyword (string fname, string key, string start="", string end="")
{
    ifstream ifs (fname);
    if (!ifs) {
        cout << "Error: cannot open file " << fname << endl;
        throw;
    }
    if (start != "")
        fstream_goto (ifs, start);
    return fstream_goto (ifs, key, 0, end);
}

vector<string> read_bracket (ifstream& ifs, string key, int skipline=0)
{
    vector<string> re;

    // Go to the key-word line
    bool good = fstream_goto (ifs, key);
    if (!good) return re;
    good = fstream_goto (ifs, "{", skipline);
    if (!good) return re;

    string line;
    while (getline(ifs, line))
    {
        auto n = line.find ("}");
        if (n != string::npos) break;

        re.push_back (line);
    }
    return re;
}

vector<string> read_bracket (string file, string key, int skipline=0)
{
    ifstream ifs (file);
    if (!ifs)
    {
        cout << "Error: read_bracket: cannot open file: " << file << endl;
        throw;
    }
    return read_bracket (ifs, key, skipline);
}

template <typename Type1, typename Type2>
vector <tuple<Type1, Type2>> read_bracket_values (string file, string key, int skipline=0)
{
    auto lines = read_bracket (file, key, skipline);
    vector<tuple<Type1, Type2>> re;
    for(auto& line : lines)
    {
        rstrip (line);
        std::istringstream iss (line);
        Type1 val1;
        Type2 val2;
        iss >> val1 >> val2;
        re.emplace_back (val1, val2);
    }
    return re;
}

template <typename Type1, typename Type2, typename Type3>
vector <tuple<Type1, Type2, Type3>> read_bracket_values (string file, string key, int skipline=0)
{
    auto lines = read_bracket (file, key, skipline);
    vector<tuple<Type1, Type2, Type3>> re;
    for(auto& line : lines)
    {
        rstrip (line);
        std::istringstream iss (line);
        Type1 val1;
        Type2 val2;
        Type3 val3;
        iss >> val1 >> val2 >> val3;
        re.emplace_back (val1, val2, val3);
    }
    return re;
}

template <typename Type1, typename Type2, typename Type3, typename Type4>
vector <tuple<Type1, Type2, Type3, Type4>> read_bracket_values (string file, string key, int skipline=0)
{
    auto lines = read_bracket (file, key, skipline);
    vector<tuple<Type1, Type2, Type3, Type4>> re;
    for(auto& line : lines)
    {
        rstrip (line);
        std::istringstream iss (line);
        Type1 val1;
        Type2 val2;
        Type3 val3;
        Type4 val4;
        iss >> val1 >> val2 >> val3 >> val4;
        re.emplace_back (val1, val2, val3, val4);
    }
    return re;
}
#endif
