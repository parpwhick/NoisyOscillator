/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Utilities.h
 * Author: dawid
 *
 * Created on February 5, 2016, 5:26 PM
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <array>
#include <vector>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Eigen/Dense"

template <typename T, std::size_t N, typename U>
std::array<T,N> operator/ (std::array<T,N> v, U scalar){
    std::array<T,N> r;
    for(size_t i = 0; i < N; i++)
        r[i] = v[i] / scalar;
    return r;
}

template <typename T, std::size_t N, typename U>
std::array<T,N> operator* (std::array<T,N> v, U scalar){
    std::array<T,N> r;
    for(size_t i = 0; i < N; i++)
        r[i] = v[i] * scalar;
    return r;
}

template <typename T, std::size_t N>
std::array<T,N> operator+ (std::array<T,N> v, std::array<T,N> u){
    std::array<T,N> r;
    for(size_t i = 0; i < N; i++)
        r[i] = v[i] + u[i];
    return r;
}

template <typename T, std::size_t N>
std::array<T,N> operator* (std::array<T,N> v, std::array<T,N> u){
    std::array<T,N> r;
    for(size_t i = 0; i < N; i++)
        r[i] = v[i] * u[i];
    return r;
}

template <typename T, std::size_t N>
std::array<T,N> operator- (std::array<T,N> v, std::array<T,N> u){
    return v + u * (-1);
}

template <typename T, std::size_t N>
std::ostream& operator<< (std::ostream& out, const std::array<T, N>& v) {
    out << '[';
    for (std::size_t i = 0; i < N; i++)
        out << v[i] << ", "; 
    out << "\b\b]";
  return out;
}

template <typename T> class windowing {
public:
    
    std::vector<T> w;
    std::vector<T> buf;
    int from;
    int N;
    
    windowing(int n){
        N = n;
        w = lanczos(n);
        buf.resize(0.0, n);
        from = 0;
    }
    ~windowing();
    
    std::vector<T> lanczos(int n) {
        //    std::cerr << "Creating window with size " << n << std::endl;
        std::vector<T> w(n);
        T div = n - 1;
        T norm = 0.0;
        for (int k = 0; k < n; k++) {
            w[k] = sinc(2 * k / div - 1);
            norm += w[k];
        }
        for (int k = 0; k < n; k++) {
            w[k] /= norm;
        }
        return w;
    }

    T window_avg() {
        T avg = 0.0;
        
        for (int i = 0; i < N; i++) {
            avg += w[i] * buf[(i + from) % N];
        }

        return avg;
    }
    
    void add_value(T val){
        buf[from] = val;
        from = from + 1 % N;
    }

};

constexpr int X = 0;
constexpr int Y = 1;
constexpr int Z = 2;

typedef std::array<double, 3> staticvec;
inline double scalar(const staticvec& a, const staticvec&b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double module(const staticvec & vector) {
    double m = scalar(vector, vector);
    return std::sqrt(m);
}

inline double normalize(staticvec& vector) {
    auto m = module(vector);
    vector= {vector[X] / m, vector[Y] / m, vector[Z] / m};
    return m;
}

template<typename T> inline T square(const T & x){
    return x*x;
}

inline double sinc(double x){
    //for x > 1e-8 (sqrt[eps]), sinc(x)==1.0
    return (std::abs(x)>1e-8) ? (std::sin(x)/x) : x;
}

template <typename T, std::size_t N>
std::array<T,N> square (std::array<T,N> v){
    std::array<T,N> r;
    for(size_t i = 0; i < N; i++)
        r[i] = square(v[i]);
    return r;
}



constexpr auto zeros = staticvec{0.0,0.0,0.0};
std::string autoFileName(std::string prefix = "output");
void print_table(std::string prefix, const Eigen::MatrixXd & table);

#endif /* UTILITIES_H */

