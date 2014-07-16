#include "align/Variant.h"
#include <string>

Variant::Variant() :
    position(std::numeric_limits<decltype(position)>::max()) {
}

const Variant Variant::none;

Variant::Variant(size_t position, ReferenceString insertionString) :
    position(position),
    insertionString(insertionString) {
}

Variant::Variant(size_t position, size_t deletionLength) :
    position(position),
    deletionLength(deletionLength) {
}

Variant::Variant(size_t position, ReferenceString insertionString, size_t deletionLength) :
    position(position),
    deletionLength(deletionLength),
    insertionString(insertionString) {
}

size_t Variant::getPosition() const {
    return position;
}

const ReferenceString &Variant::getInsertionString() const {
    return insertionString;
}

size_t Variant::getDeletionLength() const {
    return deletionLength;
}

Variant::Type Variant::getType() const {
    if (seqan::length(insertionString) == 0) {
        return DELETION;
    }
    return (deletionLength == 0) ? INSERTION : COMPLEX;
}

bool Variant::operator<(const Variant &rhs) const {
    const Variant &lhs = *this;
    auto lhsPos = lhs.getPosition();
    auto rhsPos = rhs.getPosition();
    if (lhsPos < rhsPos) {
        return true;
    } else if (lhsPos > rhsPos) {
        return false;
    }
    auto lhsDel = lhs.getDeletionLength();
    auto rhsDel = rhs.getDeletionLength();
    if (lhsDel < rhsDel) {
        return true;
    } else if (lhsDel > rhsDel) {
        return false;
    }
    auto lhsIns = seqan::length(lhs.getInsertionString());
    auto rhsIns = seqan::length(rhs.getInsertionString());
    if (lhsIns < rhsIns) {
        return true;
    } else if (lhsIns > rhsIns) {
        return false;
    }
    return &lhs < &rhs;
}
