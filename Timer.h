#ifndef __TIMERS_H_CMC__
#define __TIMERS_H_CMC__
#include <ctime>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
using namespace std;

class Timer
{
    public:
        Timer (string name="",bool start_=false)
        : _name (name)
        , _t1 (0.)
        , _t2 (0.)
        {
            if (start_) start();
        }

        // Set a new starting time.
        // The time duration between the last rec() and the new start() will not be counted.
        void start ()
        {
            clock_t t = clock();
            _t1 = t - (_t2 - _t1);
            _t2 = t;
        }
        // Record the time.
        // One can keep calling more than one rec() after a start(),
        // or making new starting point by calling start().
        void    stop    ()          { _t2 = clock(); }
        // Return the total time duration between each starts and records (in seconds).
        double  t       () const    { return double(_t2 - _t1) / CLOCKS_PER_SEC; }
        void    reset   ()          { _t1 = _t2 = 0.; start(); }
        string  name    () const    { return _name; }
        void    print   () const    { auto t = clock(); cout << _name << ": " << double(t-_t1) / CLOCKS_PER_SEC << endl; }

    private:
        string  _name;
        clock_t _t1, _t2;
};

class Timers
{
    public:
        Timers ()
        : _main ("Main")
        , _last_print_time (0.)
        , _start(false)
        {}

        void    print_duration  (double dt)     { _dt = dt*60.; } // in minutes
        void    start           ()              { _main.start(); _start = true; }
        bool    print           (bool force_print = false);
        bool    started         () const        { return _start; }
        Timer&  operator[]      (string name)
        {
            if (_timers.count(name) == 0)
                this->add (name);
            return _timers[name];
        }
        void    add             (string name)
        {
            //_timers.push_back (&timer);
            _timers[name] = Timer (name);
        }

    private:
        //vector<Timer*> _timers;
        map<string, Timer>  _timers;
        Timer               _main;
        double              _dt, _last_print_time;
        bool                _start;

        string human_time (int t);
};

bool Timers :: print (bool force_print)
{
  if (_timers.size() == 0) return false;

  _main.stop();
  double t_main = _main.t();

  if (t_main - _last_print_time > _dt || force_print) {
    cout << endl;
    cout << "    Time report\n";
    cout << "        " << _main.name() << ": " << human_time(t_main) << " (" << int(t_main) << " secs)\n";

    for(const auto& item : _timers) {
      const auto& timer = item.second;
      int t = int(timer.t());
      cout << "        " << left << setw(18) << timer.name() << ": " << t << " secs, ";
      printf("%.2f",t*100./t_main);
      cout << " %\n";
    }
    cout << endl;

    _last_print_time = t_main;
    return true;
  }
  return false;
}

string Timers :: human_time (int t) // t is in seconds
{
  int days = t / 60 / 60 / 24;
  int hours = (t / 60 / 60) % 24;
  int mins = (t / 60) % 60;
  int secs = t % 60;
  string str ("");
  if (days > 0) str += to_string(days) + " d ";
  if (hours > 0) str += to_string(hours) + " h ";
  if (mins > 0) str += to_string(mins) + "m ";
  str += to_string(secs) + "s";
  return str;
}
#endif
