#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <unordered_map>
using namespace std;
using namespace Eigen;
//изучив немного теории, решила начать сначала, написать программу для нахождения координат одной антенны, затем реализовать для трех.
//интеренет говорит, что для нелинейной системы лучше подходит или обычный GD  или SGD, а я похоже запуталась, поэтому реализую просто GD


//первые две функции не отличаются от предыдущего кода, почти.Разница в том, что гиперболу я рассчитываю как разницу между расстояниями (от А до D) и (от А до Е) например, а не как разницу расстояний двух антенн до одного и того же источника.
double evclid(Vector2d point1, Vector2d point2) {
    return (point1 - point2).norm();
}

double loss_function(Vector2d A, Vector2d D, Vector2d E, Vector2d F, const vector<double>& distance_diff) {
    double mistake = 0.0;
    mistake += pow(evclid(A, D) - evclid(A, E) - distance_diff[0], 2);
    mistake += pow(evclid(A, D) - evclid(A, F) - distance_diff[1], 2);
    mistake += pow(evclid(A, E) - evclid(A, F) - distance_diff[2], 2);
    return mistake;
}



void Gradient_spusk(Vector2d& A, Vector2d D, Vector2d E, Vector2d F, const vector<double>& distance_diff, int iterations) {

    //θ(t + 1) = θ(t)−α⋅∇f(θ(t)) реализовать надо вот это
    for (int iter = 0; iter < iterations; iter++) {
        //получается  teta  это или x  или y 
        for (int coord = 0; coord < 2; coord++) {
            // тогда выходит три градиента за один проход для  x и y. 
            double grad = 0.0;
            double normA_D = (A - D).norm();
            double normA_E = (A-E).norm();
            double normA_F = (A-F).norm();

        }



    }

}