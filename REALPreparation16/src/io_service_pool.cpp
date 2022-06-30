#include "io_service_pool.hpp"

io_service_pool::io_service_pool(const unsigned concurrency) : w(make_unique<work>(*this))
{
    reserve(concurrency);
    for (unsigned i = 0; i < concurrency; i++)
    {
        emplace_back(async(launch::async, [&]()
        {
            run();
        }));
    }
}

inline void io_service_pool::wait()
{
    this->w.reset();
    for (auto& f : *this)
    {
        f.get();
    }
}