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

template <typename T, std::size_t N>
std::ostream& operator<< (std::ostream& out, const std::array<T, N>& v) {
    out << '[';
    for (std::size_t i = 0; i < N; i++)
        out << v[i] << ", "; 
    out << "\b\b]";
  return out;
}


inline double module(vec vector) {
    double m = 0.0;
    for (auto & x : vector) {
        m += x * x;
    }
    return std::sqrt(m);
}


const int X = 0;
const int Y = 1;
const int Z = 2;

template<typename T> inline T square(const T & x){
    return x*x;
}

inline double sinc(double x){
    //for x > 1e-8 (sqrt[eps]), sinc(x)==1.0
    return (std::abs(x)>1e-8) ? (std::sin(x)/x) : x;
}

template <typename T>
std::vector<T> lanczos(int n){
//    std::cerr << "Creating window with size " << n << std::endl;
    std::vector<T> w(n);
    T div = n-1;
    T norm = 0.0;
    for (int k = 0; k < n; k++){
        w[k] = sinc(2*k/div - 1) ;
        norm += w[k];
    }
    for (int k = 0; k < n; k++){
        w[k] /= norm;
    }
    return w;    
}

#endif /* UTILITIES_H */

