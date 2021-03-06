#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "matrix.h"

// Run an activation function on each element in a matrix,
// modifies the matrix in place
// matrix m: Input to activation function
// ACTIVATION a: function to run
void activate_matrix(matrix m, ACTIVATION a)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        double sum = 0;
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            if(a == LOGISTIC){
                m.data[i][j] = 1 / (1 + exp(-x));  
            } else if (a == RELU){
                m.data[i][j] = (x > 0) ? x : 0;
            } else if (a == LRELU){
                m.data[i][j] = (x > 0) ? x : (0.1 * x);
            } else if (a == SOFTMAX){
                m.data[i][j] = exp(x); 
            }
            sum += m.data[i][j];
        }
        if (a == SOFTMAX) {
            for (j = 0; j < m.cols; ++j) {
                m.data[i][j] /= sum; 
            }
        }
    }
}

// Calculates the gradient of an activation function and multiplies it into
// the delta for a layer
// matrix m: an activated layer output
// ACTIVATION a: activation function for a layer
// matrix d: delta before activation gradient
void gradient_matrix(matrix m, ACTIVATION a, matrix d)
{
    int i, j;
    for(i = 0; i < m.rows; ++i){
        for(j = 0; j < m.cols; ++j){
            double x = m.data[i][j];
            if (a == LOGISTIC) {
                d.data[i][j] *= x * (1 - x); 
            }
            else if (a == RELU) {
                d.data[i][j] *= (x > 0) ? 1 : 0;
            }
            else if (a == LRELU) {
                d.data[i][j] *= (x > 0) ? 1 : 0.1;
            }
            else if (a == SOFTMAX) {
                d.data[i][j] *= 1; 
            }
        }
    }
}

// Forward propagate information through a layer
// layer *l: pointer to the layer
// matrix in: input to layer
// returns: matrix that is output of the layer
matrix forward_layer(layer *l, matrix in)
{

    l->in = in;  // Save the input for backpropagation

    // TODO: fix this! multiply input by weights and apply activation function.
    matrix out = matrix_mult_matrix(l->in, l->w);
    activate_matrix(out, l->activation);  

    free_matrix(l->out);// free the old output
    l->out = out;       // Save the current output for gradient calculation
    return out;
}

// Backward propagate derivatives through a layer
// layer *l: pointer to the layer
// matrix delta: partial derivative of loss w.r.t. output of layer
// returns: matrix, partial derivative of loss w.r.t. input to layer
matrix backward_layer(layer *l, matrix delta)
{
    // 1.4.1
    // delta is dL/dy
    // modify it in place to be dL/d(xw)
    gradient_matrix(l->out, l->activation, delta);

    // 1.4.2
    // then calculate dL/dw and save it in l->dw
    free_matrix(l->dw);
    matrix xt = transpose_matrix(l->in);
    matrix dw = matrix_mult_matrix(xt, delta); 
    free_matrix(xt); 
    l->dw = dw;

    
    // 1.4.3
    // finally, calculate dL/dx and return it.
    matrix wt = transpose_matrix(l->w); 
    matrix dx = matrix_mult_matrix(delta, wt);
    free_matrix(wt);

    return dx;
}

// Update the weights at layer l
// layer *l: pointer to the layer
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_layer(layer *l, double rate, double momentum, double decay)
{
    // Calculate Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    matrix temp = axpy_matrix(-decay, l->w, l->dw); 
    matrix del_w = axpy_matrix(momentum, l->v, temp); 
    free_matrix(temp); 

    // save it to l->v
    free_matrix(l->v); 
    l->v = del_w; 

    // Update l->w
    matrix new_w = axpy_matrix(rate, del_w, l->w); 
    free_matrix(l->w);
    l->w = new_w;

    // Remember to free any intermediate results to avoid memory leaks

}

// Make a new layer for our model
// int input: number of inputs to the layer
// int output: number of outputs from the layer
// ACTIVATION activation: the activation function to use
layer make_layer(int input, int output, ACTIVATION activation)
{
    layer l;
    l.in  = make_matrix(1,1);
    l.out = make_matrix(1,1);
    l.w   = random_matrix(input, output, sqrt(2./input));
    l.v   = make_matrix(input, output);
    l.dw  = make_matrix(input, output);
    l.activation = activation;
    return l;
}

// Run a model on input X
// model m: model to run
// matrix X: input to model
// returns: result matrix
matrix forward_model(model m, matrix X)
{
    int i;
    for(i = 0; i < m.n; ++i){
        X = forward_layer(m.layers + i, X);
    }
    return X;
}

// Run a model backward given gradient dL
// model m: model to run
// matrix dL: partial derivative of loss w.r.t. model output dL/dy
void backward_model(model m, matrix dL)
{
    matrix d = copy_matrix(dL);
    int i;
    for(i = m.n-1; i >= 0; --i){
        matrix prev = backward_layer(m.layers + i, d);
        free_matrix(d);
        d = prev;
    }
    free_matrix(d);
}

// Update the model weights
// model m: model to update
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_model(model m, double rate, double momentum, double decay)
{
    int i;
    for(i = 0; i < m.n; ++i){
        update_layer(m.layers + i, rate, momentum, decay);
    }
}

// Find the index of the maximum element in an array
// double *a: array
// int n: size of a, |a|
// returns: index of maximum element
int max_index(double *a, int n)
{
    if(n <= 0) return -1;
    int i;
    int max_i = 0;
    double max = a[0];
    for (i = 1; i < n; ++i) {
        if (a[i] > max){
            max = a[i];
            max_i = i;
        }
    }
    return max_i;
}

// Calculate the accuracy of a model on some data d
// model m: model to run
// data d: data to run on
// returns: accuracy, number correct / total
double accuracy_model(model m, data d)
{
    matrix p = forward_model(m, d.X);
    int i;
    int correct = 0;
    for(i = 0; i < d.y.rows; ++i){
        if(max_index(d.y.data[i], d.y.cols) == max_index(p.data[i], p.cols)) ++correct;
    }
    return (double)correct / d.y.rows;
}

// Calculate the cross-entropy loss for a set of predictions
// matrix y: the correct values
// matrix p: the predictions
// returns: average cross-entropy loss over data points, 1/n Σ(-ylog(p))
double cross_entropy_loss(matrix y, matrix p)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < y.rows; ++i){
        for(j = 0; j < y.cols; ++j){
            sum += -y.data[i][j]*log(p.data[i][j]);
        }
    }
    return sum/y.rows;
}


// Train a model on a dataset using SGD
// model m: model to train
// data d: dataset to train on
// int batch: batch size for SGD
// int iters: number of iterations of SGD to run (i.e. how many batches)
// double rate: learning rate
// double momentum: momentum
// double decay: weight decay
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay)
{
    int e;
    for(e = 0; e < iters; ++e){
        data b = random_batch(d, batch);
        matrix p = forward_model(m, b.X);
        fprintf(stderr, "%06d: Loss: %f\n", e, cross_entropy_loss(b.y, p));
        matrix dL = axpy_matrix(-1, p, b.y); // partial derivative of loss dL/dy
        backward_model(m, dL);
        update_model(m, rate/batch, momentum, decay);
        free_matrix(dL);
        free_data(b);
    }
}


// Questions 
//
// 5.2.2.1 Why might we be interested in both training accuracy and testing accuracy? What do these two numbers tell us about our current model?
// For a model, the bigger the training accuracy is, the more the model fit to the input data. 
// However, as the training accuracy grows, it migh be overfitting the traning data. 
// Thus, we will use some of the data to test the model if it only fits in some certin input. 
// Ideally, the training accruacy and the test accuracy should be similar.
//
// 5.2.2.2 Try varying the model parameter for learning rate to different powers of 10 (i.e. 10^1, 10^0, 10^-1, 10^-2, 10^-3) and training the model. What patterns do you see and how does the choice of learning rate affect both the loss during training and the final model accuracy?
// Expect for rate = 10 (which has both accuracy < 10%), all other rates are pretty accruacy (>85%). Thus all disscussion will be based on the later 4 results.
// The peak of both accuracy occur at rate = 0.1. Also, although the difference between training and testing accuracy is small, as the rate decreases, the difference between training and testing accuracy increases. 
//
// 5.2.2.3 Try varying the parameter for weight decay to different powers of 10: (10^0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5). How does weight decay affect the final model training and test accuracy?
// Before decay = 0.01, both of the accruacy increases when decay decreases. After 0.01, both accuracy remains stable. At this findal flat region, both accuracy remains similar to decay = 0.
// Different from decay = 0, where sometimes testing accuracy is bigger then training accuracy, when decay > 0, the testing accuracy is always smaller then the training accuracy.
//
// 5.2.3.1 Currently the model uses a logistic activation for the first layer. Try using a the different activation functions we programmed. How well do they perform? What's best?
// Besides softmax (which has both accuracy < 10%), all other logistic activation has a accruate result (> 90%).
// Accruacy order: RELU > LRELU > LOGISTIC > LINEAR >>>>>> SOFTMAX
//
// 5.2.3.2 Using the same activation, find the best (power of 10) learning rate for your model. What is the training accuracy and testing accuracy?
// The best rate = 0.1, leading to training accuracy = 0.96115, testing accuracy = 0.9564.
//
// 5.2.3.3 Right now the regularization parameter `decay` is set to 0. Try adding some decay to your model. What happens, does it help? Why or why not may this be?
// By setting the decay to 0.0001, the training accuracy = 0.96245 and the testing accuracy = 0.9586, which is better then setting decay = 0
// We know that neural networks are powerful models, which migh lead to overfitting. So by introducing a weight decay, we can smoothen the extreme weight comming from noise and prevent overfitting.
//
// 5.2.3.4 Modify your model so it has 3 layers instead of two. The layers should be `inputs -> 64`, `64 -> 32`, and `32 -> outputs`. Also modify your model to train for 3000 iterations instead of 1000. Look at the training and testing error for different values of decay (powers of 10, 10^-4 -> 10^0). Which is best? Why?
// By setting rate = 0.1 and decay = 0.0001 the model will be the best, with training accuracy = 0.987 and testing accuracy = 0.9751. 
// From the prevous problem, we know that adding decay can help us prevent overfitting. However, if the decay is too big, the model might be too weak. 
// By testing decay = 1, 0.1, 0.01, 0.001, 0.0001, it is seen that the smaller the decay is, the more accurate the model is. 
// This shown that the big decays will weaken the model too much, so among those five decays, we will want to use 0.0001.
//
// 5.3.2.1 How well does your network perform on the CIFAR dataset?
// I set the layer and the number of iterations as same as 5.2.3.4 and test all combination of decay = [0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1] and rate = [10, 1, 0.1, 0.01, 0.001, 0.0001]. 
// The most accuracy model uses decay = 0.0001 and rate = 0.01. The model gives us training accuracy = 0.48196, tesring accuracy = 0.4664
//



