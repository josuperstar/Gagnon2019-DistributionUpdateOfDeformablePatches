//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
template<typename T>
Vec3<T>::Vec3()
  : x(0), y(0), z(0)
{}

template<typename T>
Vec3<T>::Vec3(const T& x, const T& y, const T& z)
    : x(x), y(y), z(z)
{}

template<typename T>
Vec3<T>::Vec3(const T& xyz)
    : x(xyz), y(xyz), z(xyz)
{}

template<typename T>
Vec3<T>::Vec3(const T xyz[])
    : x(xyz[0]), y(xyz[1]), z(xyz[2])
{}

template<typename T>
Vec3<T>::Vec3(const Vec3<T>& o)
    : x(o.x), y(o.y), z(o.z)
{}

//------------------------------------------------------------------------------
// Math operators with self assignment
//------------------------------------------------------------------------------
template<typename T>
inline
Vec3<T>& Vec3<T>::operator=(const Vec3<T>& o)
{
    x = o.x;
    y = o.y;
    z = o.z;

    return *this;
}

template<typename T>
inline
Vec3<T>& Vec3<T>::operator+=(const Vec3<T> &o)
{
    x += o.x;
    y += o.y;
    z += o.z;

    return *this;
}

template<typename T>
inline
Vec3<T>& Vec3<T>::operator-=(const Vec3<T> &o)
{
    x -= o.x;
    y -= o.y;
    z -= o.z;

    return *this;
}

template<typename T>
inline
Vec3<T>& Vec3<T>::operator*=(const T& scalar)
{
    x *= scalar;
    y *= scalar;
    z *= scalar;

    return *this;
}

template<typename T>
inline
Vec3<T>& Vec3<T>::operator/=(const T& scalar)
{
    x /= scalar;
    y /= scalar;
    z /= scalar;

    return *this;
}

//------------------------------------------------------------------------------
// Math operators without self assignment (Create a temporary variable)
//------------------------------------------------------------------------------
template<typename T>
inline
Vec3<T> Vec3<T>::operator+(const Vec3<T>& o) const
{
    Vec3<T> result(*this);
    result += o;

    return result;
}

template<typename T>
inline
Vec3<T> Vec3<T>::operator-(const Vec3<T>& o) const
{
    Vec3<T> result(*this);
    result -= o;

    return result;
}

template<typename T>
inline
Vec3<T> Vec3<T>::operator*(const T& scalar) const
{
    Vec3<T> result(*this);
    result *= scalar;

    return result;
}

template<typename T>
inline
Vec3<T> Vec3<T>::operator/(const T& scalar) const
{
    Vec3<T> result(*this);
    result /= scalar;

    return result;
}

template<typename T>
inline
Vec3<T> Vec3<T>::operator-() const
{
    return Vec3<T>(-x, -y, -z);
}

template<typename T>
inline
Vec3<T> operator*(const T& scalar, const Vec3<T>& v)
{
    Vec3<T> result(v);
    result *= scalar;

    return result;
}

template<typename T>
inline
Vec3<T> operator/(const T& scalar, const Vec3<T>& v)
{
    Vec3<T> result(v);
    result /= scalar;

    return result;
}

//------------------------------------------------------------------------------
// Boolean operators
//------------------------------------------------------------------------------
template<typename T>
inline
bool Vec3<T>::operator==(const Vec3<T>& o) const
{
    return (x==o.x) && (y==o.y) && (z==o.z);
}

template<typename T>
inline
bool Vec3<T>::operator!=(const Vec3<T>& o) const
{
    return !(*this==o);
}

//------------------------------------------------------------------------------
// Vector operations
//------------------------------------------------------------------------------
template<typename T>
inline
T Vec3<T>::length2() const
{
    return x*x + y*y + z*z;
}

template<typename T>
inline
T Vec3<T>::length() const
{
    return std::sqrt(length2());
}

template<typename T>
inline
Vec3<T>& Vec3<T>::normalize()
{
    T len= length();

    (*this) /= len;

    return *this;
}

template<typename T>
inline
Vec3<T> Vec3<T>::getNormalizedVec() const
{
    Vec3<T> result(*this);
    result.Normalize();

    return result;
}

template<typename T>
inline
T Vec3<T>::dot(const Vec3<T>& o) const
{
    return x*o.x + y*o.y + z*o.z;
}

template<typename T>
inline
Vec3<T> Vec3<T>::cross(const Vec3<T>& o) const
{
    Vec3<T> result;
    result.x = y*o.z - z*o.y;
    result.y = z*o.x - x*o.z;
    result.z = x*o.y - y*o.x;

    return result;
}

template<typename T>
inline
Vec3<T>& Vec3<T>::applyCross(const T& o)
{
    T x0 = x;
    T y0 = y;

    x = y*o.z - z*o.y;
    y = z*o.x - x0-o.z;
    z = x0*o.y - y0*o.x;

    return *this;
}

template<typename T>
inline
T Vec3<T>::distance(const Vec3<T>& v1, const Vec3<T>& v2)
{
    return (v1 - v2).length();
}

//------------------------------------------------------------------------------
// Stream output
//------------------------------------------------------------------------------
template<typename T>
inline
std::ostream& operator<<(std::ostream& out, const Vec3<T> v)
{
    out << "[" << v.x << ", " << v.y << ", " << v.z << "]";
    return out;
}
