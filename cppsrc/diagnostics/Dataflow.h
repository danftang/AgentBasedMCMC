//
// Created by daniel on 24/10/2021.
//

#ifndef GLPKTEST_DATAFLOW_H
#define GLPKTEST_DATAFLOW_H
////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the dataflow pipeline concept.
// A dataflow pipeline here consists of a single producer, any number of transforms
// and a single consumer. However, any of these components can start other pipelines
// and in this way an acyclic graph can be built up (technically any acyclic graph if
// we allow the addition of a "universal producer" and a "universal consumer" like in
// minimum flow algorithms).
//
// CONCEPTS
// --------
// A producer is any object that can be called with no arguments to return a data object
// If a producer returns a const reference, the referred-to object must persist until
// the next call.
// A consumer is any object that can be called with a const reference to a data object and
// returns a boolean if the consumer is willing to accept more data.
// A transformer is any object that can be called with an rvalue reference to the consumer
// on it's downstream side and returns a consumer for its upstream side.
//
// TODO: Better to have a consumer be an object that takes a producer, and a transform
//  be an object that takes a producer and returns a producer? This makes type inference
//  much easier (since a function with no arguments must have a well defined return type)!
//
////////////////////////////////////////////////////////////////////////////////////////
namespace dataflow {
//    template<typename T> using Consumer = std::function<bool(T)>;
//
//    template<typename T> using Producer = std::function<T()>;

//    template<typename OPT>
//    static constexpr bool is_std_optional = false;
//
//    template<typename T>
//    static constexpr bool is_std_optional<std::optional<T>> = true;

//    template<typename T, typename=void>
//    static constexpr bool is_producer = false;
//    template<typename T>
//    static constexpr bool is_producer<T,std::enable_if_t<std::is_invocable_v<T>>> = is_std_optional<std::invoke_result_t<T>>;


//    template<typename DATAIN, typename DATAOUT> using Transform = std::function<Consumer<DATAIN>(Consumer<DATAOUT>)>;

    // Transform concept: any object for which OBJ(Consumer<OUT> &&) evaluates to RET such that RET(IN) evaluates to bool

//    template<typename OPT>
//    struct enable_if_std_optional { };
//
//    template<typename T>
//    struct enable_if_std_optional<std::optional<T>> {
//        typedef T type;
//    };
//    template<typename OPT>
//    using enable_if_std_optional_t = typename enable_if_std_optional<OPT>::type;


    template<typename T, typename FUNCTIONTEST=void>
    struct unary_function_traits {
    };

    template<typename IN, typename OUT>
    struct unary_function_traits<std::function<OUT(IN)>> {
        typedef void    is_valid_tag;
        typedef IN      argument_type;
        typedef OUT     result_type;
    };

    template<typename T>
    struct unary_function_traits<T, std::void_t<decltype(std::function(std::declval<T>()))>>:
        public unary_function_traits<decltype(std::function(std::declval<T>()))> { };


//    template<typename P, typename PRODUCERTEST=void>
//    struct enable_if_producer { };
//
//    template<typename P>
//    struct enable_if_producer<P, std::void_t<decltype(std::declval<P>()())>> {
//        typedef decltype(std::declval<P>()()) type;
//    };
//
//    template<typename PRODUCER>
//    using enable_if_producer_t = typename enable_if_producer<PRODUCER>::type;
//
//    template<typename T, typename TEST=void>
//    struct enable_if_consumer {
//    };
//
//    template<typename CONSUMER>
//    struct enable_if_consumer<CONSUMER, std::enable_if_t<
//            std::is_same_v<typename unary_function_traits<CONSUMER>::result_type, bool>
//            >> {
//        typedef typename unary_function_traits<CONSUMER>::argument_type type;
//    };
//
//    template<typename CONSUMER>
//    using enable_if_consumer_t = typename enable_if_consumer<CONSUMER>::type;
//
//    template<typename CONSUMER, typename TYPE>
//    using enable_if_consumer_of = std::enable_if_t<std::is_same_v<enable_if_consumer_t<CONSUMER>,TYPE>>;

//    template<typename T, typename IN, typename = void>
//    static constexpr bool is_consumer_of = false;
//    template<typename T, typename IN>
//    static constexpr bool is_consumer_of<T, IN, std::void_t<decltype(std::declval<T>()(std::declval<IN>()))>> = std::is_same_v<decltype(std::declval<T>()(std::declval<IN>())),bool>;



    // requires that TRANSFORM(DOWNSTREAMCONSUMER) is valid
    template<typename TRANSFORM, typename DOWNSTREAMCONSUMER, typename = std::void_t<typename TRANSFORM::transform_tag>>
    auto operator >>=(TRANSFORM &&transform, DOWNSTREAMCONSUMER &&downstreamConsumer) {
        return [downstreamConsumer = std::forward<DOWNSTREAMCONSUMER>(downstreamConsumer),
                transform = std::forward<TRANSFORM>(transform)](const auto &item) mutable {
            return transform(downstreamConsumer, item);
        };
    }

    // requires that CONSUMER(PRODUCER) returns bool
    template<typename PRODUCER, typename CONSUMER, typename = std::void_t<
                std::enable_if_t<
                    std::is_same_v<decltype(std::declval<CONSUMER>()(std::declval<PRODUCER>()())), bool>
                >
            >>
    void operator >>=(PRODUCER &&producer, CONSUMER &&consumer) {
        while(consumer(producer())) {};
    }

//    template<typename T>
//    class Split: public Consumer<T> {
//    public:
//        Split(std::initializer_list<Consumer<T>> consumers):
//                Consumer<T>([consumers = std::vector(consumers), isOpen = std::vector<bool>(consumers.size(),true)](T item) mutable {
//                    bool anyOpen = false;
//                    for(int i=0; i<consumers.size(); ++i) {
//                        if(isOpen[i]) {
//                            anyOpen |= (isOpen[i] = isOpen[i] & consumers[i](item));
//                        }
//                    }
//                    return anyOpen;
//                }) {}
//    };

    class Transform {
    public:
        typedef void transform_tag;
    };

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


    class Take: public Transform {
    public:
        int countDown;
        Take(int n) : countDown(n) { }

        template<typename CONSUMER, typename T>
        bool operator()(CONSUMER &downstreamConsumer, const T &item) {
            bool downstreamOpen = downstreamConsumer(item);
            --countDown;
            return (downstreamOpen && countDown != 0);
        }
    };

    template<typename FUNC, typename = typename unary_function_traits<FUNC>::is_valid_tag>
    class Map: public Transform {
    public:
        FUNC mapFunc;

        typedef typename unary_function_traits<FUNC>::argument_type upstream_type;

        Map(FUNC mapFunction): mapFunc(mapFunction) { }

        // requires that
        template<typename CONSUMER>
        bool operator()(CONSUMER &downstreamConsumer, const upstream_type &item) {
            return downstreamConsumer(mapFunc(item));
        }
    };


    class Drop: public Transform {
    public:
        int n;
        Drop(int n): n(n) {}

        template<typename CONSUMER, typename T>
        bool operator()(CONSUMER &downstreamConsumer, const T &item) {
            if(n>0) {
                --n;
                return true;
            }
            return downstreamConsumer(item);
        }
    };

    // passes 1 then drops n-1 then repeats
    class Thin: public Transform {
    public:
        const int n;
        int count=0;

        Thin(int n): n(n) {}

        template<typename CONSUMER, typename T>
        bool operator()(CONSUMER &downstreamConsumer, const T &item) {
            if(count%n == 0) {
                count = 1;
                return downstreamConsumer(item);
            }
            ++count;
            return true;
        }
    };

    // TODO: can we get type inference here?
    //  Would it require a separation of connection and data processing steps, or could we use template templates?
    template<typename DATA>
    class CollectThenEmit: public Transform {
    public:
        const int n;
        std::vector<DATA> data;
        CollectThenEmit(int n): n(n) {
            data.reserve(n);
        }

        template<typename CONSUMER>
        bool operator()(CONSUMER &downstreamConsumer, const DATA &item) {
            data.push_back(item);
            if(data.size() == n) {
                bool wantMore = downstreamConsumer(data);
                data.clear();
                return wantMore;
            }
            return true;
        }

    };

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

//    template<typename FIRSTCONSUMER,typename LASTCONSUMER>
//    class SwitchOnClose {
//    public:
//        FIRSTCONSUMER consumer1;
//        LASTCONSUMER  consumer2;
//        bool consumer1IsOpen;
//
//        SwitchOnClose(FIRSTCONSUMER beforeSwitch, LASTCONSUMER afterSwitch):
//                consumer1(std::move(beforeSwitch)),
//                consumer2(std::move(afterSwitch)),
//                consumer1IsOpen(true) {}
//
//        template<typename T>
//        bool operator()(const T &item) {
//            if (consumer1IsOpen) {
//                consumer1IsOpen = consumer1(item);
//                return true;
//            }
//            return consumer2(item);
//        }
//    };


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


    template<typename DATA>
    class Sum {
    public:
        DATA &sum;

        Sum(DATA &result): sum(result) {}

        bool operator()(const DATA &item) {
            sum += item;
            return true;
        }
    };


    template<typename T>
    auto save(std::vector<T> &vectorLog) {
        return [&vectorLog](const T &item) {
            vectorLog.push_back(item);
            return true;
        };
    }

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
