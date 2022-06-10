#ifndef IO_SERVICE_POOL_HPP
#define IO_SERVICE_POOL_HPP
#include <iostream>
#include <future>
#include <boost/asio/io_service.hpp>
#include <vector>

using namespace std;
using namespace boost::asio;

class io_service_pool : public io_service, public vector<future<void>>
{
public:
    explicit io_service_pool(const unsigned);
    inline void wait();
private:
    unique_ptr<work> w;
};


#endif