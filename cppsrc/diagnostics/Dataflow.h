//
// Created by daniel on 24/10/2021.
//

#ifndef GLPKTEST_DATAFLOW_H
#define GLPKTEST_DATAFLOW_H

namespace dataflow {
    template<typename T> using Consumer = std::function<bool(T)>;

    template<typename T> using Producer = std::function<std::optional<T>()>;

//    template<typename DATAIN, typename DATAOUT> using Transform = std::function<Consumer<DATAIN>(Consumer<DATAOUT>)>;

    // Transform concept: any object for which OBJ(Consumer<OUT> &&) evaluates to RET such that RET(IN) evaluates to bool

    template<typename TRANSFORMCONCEPT, typename OUT>
    auto operator >>=(TRANSFORMCONCEPT &&t, dataflow::Consumer<OUT> downStreamConsumer) {
        return t(std::move(downStreamConsumer));
    }

//    template<typename TRANSFORMCONCEPT, typename OUT>
//    auto operator >>=(TRANSFORMCONCEPT &&t, dataflow::Consumer<OUT> &downStreamConsumer) {
//        Consumer<OUT> moveableDownstream = [&downStreamConsumer](OUT item) { return downStreamConsumer(item);};
//        return t(std::move(moveableDownstream));
//    }

    template<typename T>
    void operator >>=(dataflow::Producer<T> &producer, dataflow::Consumer<T> consumer) {
        std::optional<T> nextItem = producer();
        bool channelIsOpen = nextItem.has_value();
        while(channelIsOpen) {
            channelIsOpen = consumer(nextItem.value());
            nextItem = producer();
            channelIsOpen &= nextItem.has_value();
        }
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

    template<typename SIG>
    class Map {
    public:
        std::function<SIG> mapFunc;
        typedef typename std::function<SIG>::result_type out_type;
        typedef typename std::function<SIG>::argument_type in_type;

        Map(std::function<SIG> mapFunction): mapFunc(mapFunction) { }

        Consumer<in_type> operator()(Consumer<out_type> downstreamConsumer) {
            return [downstreamConsumer = std::move(downstreamConsumer), mapFunc=this->mapFunc](in_type item) {
                return downstreamConsumer(mapFunc(item));
            };
        }
    };

    template<typename T>
    Consumer<T> switchAfter(int n, Consumer<T> beforeSwitch, Consumer<T> afterSwitch) {
        return [n, beforeSwitch=std::move(beforeSwitch), afterSwitch=std::move(afterSwitch)](T item) mutable {
            if(n > 0) {
                --n;
                return beforeSwitch(item);
            }
            return afterSwitch(item);
        };
    }



//    template<typename... ARGS>
//    auto consumer(ARGS &&... args) {
//        return (std::move(args) >>= ...);
//    }

//    template<typename IN, typename OUT>
//    Transform<IN,OUT> map(std::function<OUT(IN)> t) {
//        return [outert=std::move(t)](Consumer<OUT> downStreamConsumer) {
//            return [innert=std::move(outert), downStreamConsumer](IN inItem) { return downStreamConsumer(innert(inItem)); };
//        };
//    }


//    template<typename LEFT, typename FUNC>
//    auto operator >> (Producer<LEFT> &&p, FUNC map) -> Producer<decltype(map(p().value()))> {
//        typedef decltype(map(p().value())) RETURN;
//        return [&p,map]() {
//            const std::optional<LEFT> &nextItem = p();
//            return nextItem.has_value()?std::optional(map(nextItem.value())):std::optional<RETURN>();
//        };
//    }
//
//
//
//    template<typename LEFT, typename RIGHT>
//    Producer<RIGHT> operator >> (Producer<LEFT> &&p, std::function<std::optional<RIGHT>(std::optional<LEFT>)> map) {
//        return [&p,map]() {
//            return map(p());
//        };
//    }



//    template<typename LEFT, typename MID, typename RIGHT>
//    Transform<LEFT,RIGHT> operator >>(Transform<LEFT,MID> leftTransform, Transform<MID,RIGHT> rightTransform) {
//        return [leftTransform = std::move(leftTransform),rightTransform=std::move(rightTransform)](Consumer<RIGHT> &&c) {
//            return leftTransform(rightTransform(c));
//        };
//    }
//
//    template<typename LEFT, typename MID, typename FUNC>
//    auto operator >>(Transform<LEFT,MID> leftTransform, FUNC map) -> Transform<LEFT,decltype(map(*static_cast<MID*>(nullptr)))> {
//        typedef decltype(map(*static_cast<MID*>(nullptr))) RIGHT;
//        return [leftTransform=std::move(leftTransform), map](Consumer<RIGHT> downStreamConsumer) {
//            return leftTransform([map, downStreamConsumer](MID inItem) { return downStreamConsumer(map(inItem)); });
//        };
//    }
//

//        template<typename T>
//        operator Transform<T,T>() {
//            return [n=countDown](Consumer<T> c) {
//                return [n, c](T item) mutable {
//                    bool downstreamOpen = c(item);
//                    --n;
//                    return (downstreamOpen && n != 0);
//                };
//            };
//        }


//    template<typename T>
//    Producer<T> operator >>(Producer<T> &p, Take &&t) {
//        return [&p,&t]() mutable {
//            if(t.countDown == 0) return std::optional<T>();
//            --t.countDown;
//            return p();
//        };
//    }
//
//    template<typename T>
//    Consumer<T> operator<<(Consumer<T> &c, Take &&t) {
//        return [&c,&t](T item) {
//            bool downstreamOpen = c(item);
//            --t.countDown;
//            return (downstreamOpen && t.countDown != 0);
//        };
//    }

//    template<typename T>
//    Transform<T,T> takeTransform(int n) {
//        return [n](Consumer<T> c) {
//            return [n, c](T item) mutable {
//                bool downstreamOpen = c(item);
//                --n;
//                return (downstreamOpen && n != 0);
//            };
//        };
//    }

};

//template<typename T>
//void operator |(dataflow::Producer<T> producer, dataflow::Consumer<T> consumer) {
//    std::optional<T> nextItem = producer();
//    bool channelIsOpen = nextItem.has_value();
//    while(channelIsOpen) {
//        channelIsOpen = consumer(nextItem.value());
//        nextItem = producer();
//        channelIsOpen &= nextItem.has_value();
//    }
//}

//template<typename PRODUCED, typename CHAIN>
//Chain<CHAIN> operator >>(dataflow::Producer<PRODUCED> &producer, dataflow::Transform<PRODUCED,CHAIN> transform) {
//    return Chain(producer, transform);
//}

//template<typename IN, typename OUT>
//dataflow::Consumer<IN> operator >>(dataflow::Transform<IN,OUT> t, dataflow::Consumer<OUT> &&downStreamConsumer) {
//    return t(std::move(downStreamConsumer));
//}




#endif //GLPKTEST_DATAFLOW_H