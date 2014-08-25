#ifndef PG583_ALIGN_VARIANT_H
#define PG583_ALIGN_VARIANT_H

#include "types.h"

#include <stddef.h>

class Variant {
    Variant(); // private, only used for static instance "none"
public:
    static const Variant none;

    enum Type {
        INSERTION, DELETION, COMPLEX
    };

    Variant(size_t position, size_t deletionLength);
    Variant(size_t position, const ReferenceString &insertionString);
    Variant(size_t position, const ReferenceString &insertionString, size_t deletionLength);

    // move constructor and assignment operator need for std::sort
    Variant(Variant &&) = default;
    Variant& operator=(Variant &&var) = default;

    Variant(const Variant &) = default;
    Variant& operator=(const Variant &var) = default;
    ~Variant() = default;

    size_t getPosition() const;
    const ReferenceString& getInsertionString() const;
    size_t getDeletionLength() const;

    Type getType() const;

    // needed for std::set and std::sort
    bool operator<(const Variant &rhs) const;
private:
    size_t position;
    size_t deletionLength = 0;
    ReferenceString insertionString = "";
};

// needed for containers in Alignment.h and BacktracingMatrix.h
template<class T>
bool operator<(const std::reference_wrapper<T> &lhs, const std::reference_wrapper<T> &rhs) {
    return (T &) lhs < (T &) rhs;
}

#endif // PG583_ALIGN_VARIANT_H
