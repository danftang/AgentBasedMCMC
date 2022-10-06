//
// Created by daniel on 24/10/2021.
//

#ifndef GLPKTEST_DATAFLOW_H
#define GLPKTEST_DATAFLOW_H
////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the dataflow pipeline concept.
//
// CONCEPTS
// --------
// A PRODUCER is any object that can be called with no arguments to return a data object
// If a producer returns a const reference, the referred-to object must persist until
// the next call.
//
// A CONSUMER is any object that can be called with a const reference to a data object and
// returns a boolean if the consumer is willing to accept more data.
//
// A TRANSFORMER is any object that can be called with a const referene to an upstream data
// item, and returns a data item that is to be fed downstream. Any object that can
// be used to construct a std::function<OUt(IN)> object can automatically be used as
// a transformer. A member function of a downstream data item can also be used as a transformer.
// Transformers with indeterminate upstream/downstream type can be
// created with a class that supplies a templated operator () and an overriden >>=
// operator.
//
// Elements are connected together using the '=>' operator. A dataflow PIPELINE consists of
// a PRODUCER >>= CONSUMER pair, where a TRANSFORMER >>= CONSUMER pair is also a CONSUMER.
// By allowing components to be constructed from other CONSUMERS a tree-like
// PIPELINE can be constructed.
//
// So for example, the code:
//
// producer => transform => consumer;
//
// would describe the typical pipeline.
//
////////////////////////////////////////////////////////////////////////////////////////
#include <future>
#include "../gnuplot-iostream/gnuplot-iostream.h"

namespace dataflow {

    // To maintain c++20 compatability...
    template<class T> struct unary_function_arg { };

    template<class IN, class OUT>
    struct unary_function_arg<std::function<OUT(IN)>> {
        typedef IN type;
    };

    template<class CONSUMER, class ITEM>
    using enable_if_is_consumer_of = std::enable_if_t<std::is_invocable_r_v<bool,CONSUMER,ITEM>>;

    // Transform by object that can be constructor of a std::function<OUT(IN)>
    template<class TRANSFORM, class CONSUMER,
            class FUNC = decltype(std::function(std::declval<TRANSFORM>())),
            class ITEM = typename unary_function_arg<FUNC>::type,
            typename = enable_if_is_consumer_of<CONSUMER,typename FUNC::result_type>>
    auto operator >>=(TRANSFORM &&transform, CONSUMER &&consumer) {
        return [trans = std::forward<TRANSFORM>(transform), consumer = std::forward<CONSUMER>(consumer)](const ITEM &item) mutable -> bool {
            return consumer(trans(item));
        };
    }

    // Transform by data-type member function
    template<class RETURN, class ITEM, class CONSUMER,
            typename = enable_if_is_consumer_of<CONSUMER,RETURN>>
    auto operator >>=(RETURN(ITEM::*ptr)() const, CONSUMER &&consumer) {
        return [ptr, consumer = std::forward<CONSUMER>(consumer)](const ITEM &item) mutable {
            return consumer((item.*ptr)());
        };
    }

    // execute a dataflow diagram
    template<typename PRODUCER, typename CONSUMER,
            typename = enable_if_is_consumer_of<CONSUMER,std::result_of_t<PRODUCER()>>>
    inline void operator >>=(PRODUCER &&producer, CONSUMER &&consumer) {
        while(consumer(producer())) {};
    }


    // explicit exec
    template<typename PRODUCER, typename CONSUMER>
    inline void exec(PRODUCER &&producer, CONSUMER &&consumer) {
        while(consumer(producer())) {};
    }

    // execute a dataflow diagram in a separate thread
    template<typename PRODUCER, typename CONSUMER>
    inline auto exec_async(PRODUCER &&producer, CONSUMER &&consumer) {
        return std::async([producer = std::forward<PRODUCER>(producer), consumer = std::forward<CONSUMER>(consumer)]() {while(consumer(producer())) {};});
    }



    // Split is a CONSUMER that splits a datastream over any number of consumers,
    // sending a copy of each data item to each consumer until all consumers have closed
    // the connection
    // e.g. producer => Split(consumer1, consumer 2)
    template<typename... CONSUMERS>
    class Split {
    public:
        std::tuple<CONSUMERS...> consumers;
        std::vector<bool> isOpen;

        Split(CONSUMERS... consumers):
        consumers(std::move(consumers)...),
        isOpen(sizeof...(consumers), true) {}

        template<typename T>
        bool operator()(const T &item) {
            feedConsumers(item, std::index_sequence_for<CONSUMERS...>());
            return std::find(isOpen.begin(),isOpen.end(),true) != isOpen.end();
        }

    private:
        template<typename T, size_t...INDEXES>
        void feedConsumers(const T &item, std::index_sequence<INDEXES...> idx) {
            ((isOpen[INDEXES]?(isOpen[INDEXES] = std::get<INDEXES>(consumers)(item)):false),...);
        }
    };


    // Take is a TRANSFORM that takes 'n' data items and then closes the connection
    class Take {
    public:
        int countDown;
        Take(int n) : countDown(n) { }

        template<typename CONSUMER>
        auto operator >>=(CONSUMER &&downstreamConsumer) {
            return [countDown = this->countDown, downstreamConsumer](const auto &item) mutable -> bool {
                bool downstreamOpen = downstreamConsumer(item);
                --countDown;
                return (downstreamOpen && countDown != 0);
            };
        }
    };



    // Drop is a TRANSFORM that drops the first 'n' data items then passes
    // data downstream
    class Drop {
    public:
        int n;
        Drop(int n): n(n) {}

        template<typename CONSUMER>
        inline auto operator >>=(CONSUMER &&downstreamConsumer) {
            return [n = this->n, downstreamConsumer](const auto &item) mutable -> bool {
                if (n > 0) {
                    --n;
                    return true;
                }
                return downstreamConsumer(item);
            };
        }
    };

    // passes 1 then drops n-1 then repeats
    class Thin {
    public:
        const int interval;

        Thin(int interval): interval(interval) {}

        template<typename CONSUMER>
        inline auto operator >>=(CONSUMER &&downstreamConsumer) {
            return [n = interval, count = 0, downstreamConsumer](const auto &item) mutable -> bool {
                if (count % n == 0) {
                    count = 1;
                    return downstreamConsumer(item);
                }
                ++count;
                return true;
            };
        }
    };



    // Repeatedly collects 'n' data items, packages them into a vector and sends them downstream
    class CollectThenEmit {
    public:
        const int n;
        CollectThenEmit(int n): n(n) {
        }

        template<typename CONSUMER>
        auto operator >>=(CONSUMER &downstreamConsumer) {
            return [n = n, downstreamConsumer](const auto &item) mutable -> bool {
                static thread_local std::vector<decltype(std::remove_const(std::remove_reference(item)))> data; // now we know item type
                if(data.size() == 0) data.reserve(n);
                data.push_back(item);
                if (data.size() == n) {
                    bool wantMore = downstreamConsumer(data);
                    data.clear();
                    return wantMore;
                }
                return true;
            };
        }

    };

    // Sends the first 'n' data items to one consumer, then the rest to another
    template<typename FIRSTCONSUMER,typename LASTCONSUMER>
    class SwitchAfter {
    public:
        FIRSTCONSUMER consumer1;
        LASTCONSUMER  consumer2;
        int switchCountdown;
        bool consumer1IsOpen;

        SwitchAfter(int nSwitchAfter, FIRSTCONSUMER beforeSwitch, LASTCONSUMER afterSwitch):
        consumer1(std::move(beforeSwitch)),
        consumer2(std::move(afterSwitch)),
        switchCountdown(nSwitchAfter),
        consumer1IsOpen(true) {}

        template<typename T>
        bool operator()(const T &item) {
            if (switchCountdown > 0) {
                --switchCountdown;
                if(consumer1IsOpen) consumer1IsOpen = consumer1(item);
                return true;
            }
            return consumer2(item);
        }
    };


    // Sends data items to the first in an ordered list of consumers
    // until that consumer closes the connection, then sends data to
    // the next consumer in the list...etc...
    template<typename... CONSUMERS>
    class SwitchOnClose {
    public:
        std::tuple<CONSUMERS...> consumers;
        int activeConsumer;

        SwitchOnClose(CONSUMERS... consumers):
                consumers(std::move(consumers)...),
                activeConsumer(0) {
        }

        template<typename T>
        bool operator()(const T &item) {
            feedConsumers(item, std::index_sequence_for<CONSUMERS...>());
            return activeConsumer < sizeof...(CONSUMERS);
        }

    private:
        template<typename T, size_t...INDEXES>
        void feedConsumers(const T &item, std::index_sequence<INDEXES...> idx) {
            bool isStillOpen = false;
            (((activeConsumer == INDEXES)?(isStillOpen = std::get<INDEXES>(consumers)(item)):false),...);
            activeConsumer += 1 - isStillOpen;
        }
    };


    // Sums data elements to a supplied accumulator
    template<typename SUM>
    class Sum {
    public:
        int n;
        SUM &sum;

        Sum(int nItems, SUM &result): sum(result), n(nItems) {}

        template<class ITEM>
        bool operator()(const ITEM &item) {
            sum += item;
            --n;
            return n > 0;
        }
    };


    class ToOstream {
    public:
        std::ostream &out;
        std::string   delimiter;

        ToOstream(std::ostream &stream, std::string Delimiter = "\n"): out(stream), delimiter(Delimiter) {}

        template<class DATA>
        bool operator()(const DATA &item) {
            out << item << delimiter;
            return true;
        }

    };

    // saves data to a vector
    template<typename T>
    auto save(std::vector<T> &vectorLog) {
        return [&vectorLog](const T &item) {
            vectorLog.push_back(item);
            return true;
        };
    }


    // saves data to a valarray, closes connection when full
    template<typename T>
    auto save(std::valarray<T> &log) {
        return [&log,i=0](const T &item) mutable {
            assert(i < log.size());
            log[i++] = item;
            return i<log.size();
        };
    }

    template<typename T=double>
    class Plot1D {
    public:
        std::vector<T> data;
        std::string title;

        Plot1D(std::string title = "dataflow Plot1D"): title(title) { }

        ~Plot1D() {
            Gnuplot gp;
            gp << "plot '-' title '" << title << "'with lines\n";
            gp.send1d(data);
        }

        bool operator()(const T &item) {
            data.push_back(item);
            return true;
        }
    };

    template<typename T=double>
    auto plot1DAfter(int nSamples, std::string title = "dataflow plot") {
        return [dataLog = std::vector<double>(nSamples), title, index = 0](const T &item) mutable {
            if(index < dataLog.size()) {
                dataLog[index] = item;
            }
            if(++index == dataLog.size()) {
                Gnuplot gp;
                gp << "plot '-' title '" << title << "'with lines\n";
                gp.send1d(dataLog);
                return false;
            }
            return true;
        };
    }


};

#endif //GLPKTEST_DATAFLOW_H
