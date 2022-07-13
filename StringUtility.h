#ifndef __STRINGUTILITY_H_CMC__
#define __STRINGUTILITY_H_CMC__
#include <string>
using namespace std;

inline bool hasStr (string str, string key)
{
    size_t pos = str.find (key);
    return !(pos == std::string::npos);
}

inline void eraseSubStr (string& mainStr, const string& toErase)
{
	// Search for the substring in string
	size_t pos = mainStr.find(toErase); 
	if (pos != string::npos)
	{
		// If found then erase it from string
		mainStr.erase(pos, toErase.length());
	}
}

// trim from left
inline void lstrip (std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
}

// trim from right
inline void rstrip (std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
}

// trim from left & right
inline void strip (std::string& s, const char* t = " \t\n\r\f\v")
{
    rstrip(s, t);
    lstrip(s, t);
}

template <typename T>
T convert_type (const string& str)
{
    if constexpr (is_same_v <T, string>)
    {
        return str;
    }
    else
    {
        T x;
        stringstream ss;
        ss << str;
        ss >> x;
        return x;
    }
}

template <typename T=string>
inline vector<T> splitStr (string str, const string& delimiter)
{
    vector<T> re;
	size_t pos = str.find (delimiter);
    while (pos != string::npos)
    {
        T x = convert_type<T> (str.substr (0, pos));
        re.push_back (x);
        str = str.substr (pos + delimiter.size());
        pos = str.find (delimiter);
    }
    T x = convert_type<T> (str);
    re.push_back (x);
    return re;
}

template <typename T>
vector<T> split_str (const string& str, int skip=0)
{
    string str2 = str;
    rstrip (str2);
    vector<T> re;
    std::istringstream iss (str2);

    string drop;
    for(int i = 0; i < skip; i++)
        iss >> drop;

    T a;
    while (iss.good()) {
        iss >> a;
        re.push_back (a);
    }
    return re;
}

string raw_string (const string& str)
{
    // s is our escaped output string
    std::string s = "";
    // loop through all characters
    for(char c : str)
    {
        // check if a given character is printable
        // the cast is necessary to avoid undefined behaviour
        if(isprint((unsigned char)c))
            s += c;
        else
        {
            std::stringstream stream;
            // if the character is not printable
            // we'll convert it to a hex string using a stringstream
            // note that since char is signed we have to cast it to unsigned first
            stream << std::hex << (unsigned int)(unsigned char)(c);
            std::string code = stream.str();
            s += std::string("\\x")+(code.size()<2?"0":"")+code;
            // alternatively for URL encodings:
            //s += std::string("%")+(code.size()<2?"0":"")+code;
        }
    }
    return s;
}
#endif
