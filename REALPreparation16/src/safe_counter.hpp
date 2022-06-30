#pragma once
#include <iostream>
#include <mutex>
#include <condition_variable>

using namespace std;

template<typename T>
class safe_counter
{
public:
    inline void init(const T);
    inline void increment();
    inline void wait();
private:
    mutex m;
    condition_variable cv;
    T n;
    T i;
};