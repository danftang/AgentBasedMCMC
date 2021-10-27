//
// Created by daniel on 24/10/2021.
//

#ifndef GLPKTEST_DATAFLOW_H
#define GLPKTEST_DATAFLOW_H

namespace dataflow {
    template<typename T> using Consumer = std::function<bool(T)>;

    template<typename T> using Producer = std::function<T()>;

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


    template<typename P, typename PRODUCERTEST=void>
    struct enable_if_producer { };

    template<typename P>
    struct enable_if_producer<P, std::void_t<decltype(std::declval<P>()())>> {
        typedef decltype(std::declval<P>()()) type;
    };

    template<typename PRODUCER>
    using enable_if_producer_t = typename enable_if_producer<PRODUCER>::type;

    template<typename T, typename TEST=void>
    struct enable_if_consumer {
    };

    template<typename CONSUMER>
    struct enable_if_consumer<CONSUMER, std::enable_if_t<
            std::is_same_v<typename unary_function_traits<CONSUMER>::result_type, bool>
            >> {
        typedef typename unary_function_traits<CONSUMER>::argument_type type;
    };

    template<typename CONSUMER>
    using enable_if_consumer_t = typename enable_if_consumer<CONSUMER>::type;

    template<typename CONSUMER, typename TYPE>
    using enable_if_consumer_of = std::enable_if_t<std::is_same_v<enable_if_consumer_t<CONSUMER>,TYPE>>;

//    template<typename T, typename IN, typename = void>
//    static constexpr bool is_consumer_of = false;
//    template<typename T, typename IN>
//    static constexpr bool is_consumer_of<T, IN, std::void_t<decltype(std::declval<T>()(std::declval<IN>()))>> = std::is_same_v<decltype(std::declval<T>()(std::declval<IN>())),bool>;




    // requires that DOWNSTREAMCONSUMER is a consumer and TRANSFORM(DOWNSTREAMCONSUMER) is a consumer
    template<typename TRANSFORM, typename DOWNSTREAMCONSUMER, typename CONDITIONS = std::void_t<
            enable_if_consumer_t<DOWNSTREAMCONSUMER>,
            enable_if_consumer_t<decltype(std::declval<TRANSFORM>()(std::declval<DOWNSTREAMCONSUMER>()))>
            >>
    auto operator >>=(TRANSFORM &&t, DOWNSTREAMCONSUMER downStreamConsumer) {
        return t(std::move(downStreamConsumer));
    }

    // requires that PRODUCER is a producer of T and CONSUMER is a consumer of T
    template<typename PRODUCER, typename CONSUMER, typename = std::void_t<
                std::enable_if_t<
                    std::is_same_v<decltype(std::declval<CONSUMER>()(std::declval<PRODUCER>()())), bool>
                >
            >>
    void operator >>=(PRODUCER &producer, CONSUMER consumer) {
        while(consumer(producer())) {};
    }

    template<typename T>
    class Split: public Consumer<T> {
    public:
        Split(std::initializer_list<Consumer<T>> consumers):
                Consumer<T>([consumers = std::vector(consumers), isOpen = std::vector<bool>(consumers.size(),true)](T item) mutable {
                    bool anyOpen = false;
                    for(int i=0; i<consumers.size(); ++i) {
                        if(isOpen[i]) {
                            anyOpen |= (isOpen[i] = isOpen[i] & consumers[i](item));
                        }
                    }
                    return anyOpen;
                }) {}
    };


    class Take {
    public:
        const int countDown;
        Take(int n) : countDown(n) { }

        template<typename T>
        Consumer<T> operator()(Consumer<T> downstreamConsumer) {
            return [c = std::move(downstreamConsumer), n=countDown](T item) mutable {
                bool downstreamOpen = c(item);
                --n;
                return (downstreamOpen && n != 0);
            };
        }
    };

    template<typename FUNC, typename = typename unary_function_traits<FUNC>::is_valid_tag>
    class Map {
    public:
        FUNC mapFunc;

        typedef typename unary_function_traits<FUNC>::argument_type upstream_type;
        typedef typename unary_function_traits<FUNC>::result_type   downstream_type;

        Map(FUNC mapFunction): mapFunc(mapFunction) { }

        template<typename CONSUMER>
        Consumer<upstream_type> operator()(CONSUMER downstreamConsumer) {
            return [downstreamConsumer = std::move(downstreamConsumer), mapFunc=this->mapFunc](upstream_type item) {
                return downstreamConsumer(mapFunc(item));
            };
        }
    };

    class Drop {
    public:
        const int n;
        Drop(int n): n(n) {}

        template<typename T>
        Consumer<T> operator()(Consumer<T> downstreamConsumer) {
            return [n=this->n, downstreamConsumer = std::move(downstreamConsumer)](T item) mutable {
                if(n>0) {
                    --n;
                    return true;
                }
                return downstreamConsumer(item);
            };
        }
    };

    template<typename T>
    class SwitchAfter: public Consumer<T> {
    public:
        SwitchAfter(int n, Consumer<T> beforeSwitch, Consumer<T> afterSwitch):
        Consumer<T>([n, beforeSwitch = std::move(beforeSwitch), afterSwitch = std::move(afterSwitch)](T item) mutable {
            if (n > 0) {
                --n;
                return beforeSwitch(item);
            }
            return afterSwitch(item);
        }) { }
    };

    template<typename T>
    Consumer<const T &> pushBack(std::vector<T> &vectorLog) {
        return [&vectorLog](const T &item) {
            vectorLog.push_back(item);
            return true;
        };
    }
};

#endif //GLPKTEST_DATAFLOW_H
