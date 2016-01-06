#ifndef OPTIMIZATIONUTILITY_HPP
#define OPTIMIZATIONUTILITY_HPP

#include <iostream>
#include <cmath>

template <typename T>
struct Vector3D {

    Vector3D(T _x, T _y, T _z) : x(_x), y(_y), z(_z) { }

    Vector3D(const Vector3D<T>& other) {

        x = other.x;
        y = other.y;
        z = other.z;
    }

    Vector3D<T>& operator= (const Vector3D<T>& rhs) {

        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    Vector3D<T>& operator+=(const Vector3D<T>& rhs) {

        x = x + rhs.x;
        y = y + rhs.y;
        z = z + rhs.z;
        return *this;
    }

    Vector3D<T>& operator-=(const Vector3D<T>& rhs) {

        x = x - rhs.x;
        y = y - rhs.y;
        z = z - rhs.z;
        return *this;
    }

    Vector3D<T>& operator/=(T rhs) {

        x = x / rhs;
        y = y / rhs;
        z = z / rhs;
        return *this;
    }

    Vector3D<T>& operator*=(T rhs) {

        x = x * rhs;
        y = y * rhs;
        z = z * rhs;
        return *this;
    }

    T dot(const Vector3D<T>& vec) const {

        return x*vec.x + y*vec.y + z*vec.z;
    }

    Vector3D<T> cross(const Vector3D<T>& vec) const {

        return Vector3D<T>(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x);
    }

    T norm() const {

        return sqrt(dot(*this));
    }

    T squared_norm() const {

        return dot(*this);
    }

    void normalize() {

        *this /= T(norm());
    }

    Vector3D<T> normalized() const {

        return Vector3D<T>(*this) /= norm();
    }

    bool is_zero() const {

        return ( norm < 0.0000000001 );
    }

    T x, y, z;
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const Vector3D<T>& vec) {
    out << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return out;
}

template <typename T>
Vector3D<T> operator+(const Vector3D<T>& left, const Vector3D<T>& right) {
    return Vector3D<T>(left.x + right.x, left.y + right.y, left.z + right.z);
}

template <typename T>
Vector3D<T> operator-(const Vector3D<T>& vec) {
    return T(-1)*vec;
}

template <typename T>
Vector3D<T> operator-(const Vector3D<T>& left, const Vector3D<T>& right) {
    return Vector3D<T>(left.x - right.x, left.y - right.y, left.z - right.z);
}

template <typename T>
Vector3D<T> operator/(const Vector3D<T>& left, T right) {
    return Vector3D<T>(left.x / right, left.y / right, left.z / right);
}

template <typename T>
Vector3D<T> operator*(const Vector3D<T>& left, T right) {
    return Vector3D<T>(left.x * right, left.y * right, left.z * right);
}

template <typename T>
Vector3D<T> operator*(T left, const Vector3D<T>& right) {
    return Vector3D<T>(right.x * left, right.y * left, right.z * left);
}

#endif // OPTIMIZATIONUTILITY_HPP
