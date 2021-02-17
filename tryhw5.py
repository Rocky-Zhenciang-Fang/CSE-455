from uwimg import *

def softmax_model(inputs, outputs):
    l = [make_layer(inputs, outputs, SOFTMAX)]
    return make_model(l)

def neural_net(inputs, outputs):
    print(inputs)
    l = [   make_layer(inputs, 64, RELU), 
            make_layer(64, 32, RELU),
            make_layer(32, outputs, SOFTMAX)]
    return make_model(l)

print("loading data...")
train = load_classification_data(b"cifar.train", b"cifar/labels.txt", 1)
test  = load_classification_data(b"cifar.test", b"cifar/labels.txt", 1)
print("done")
print

decays = [0.1, 0.001, 0.0001, 0.00001]
rates = [10, 1, 0.1, 0.01, 0.001, 0.0001]
results = []

for d in decays:
    for r in rates:
        print("training model...")
        batch = 128
        iters = 3000
        rate = r
        momentum = .9
        decay = d

        m = neural_net(train.X.cols, train.y.cols)
        train_model(m, train, batch, iters, rate, momentum, decay)
        print("done")
        print

        print("evaluating model...")
        print("decay = " + str(decay) + ", rate = " + str(rate))
        training = accuracy_model(m, train)
        testing = accuracy_model(m, test)
        print("training accuracy: %f", training)
        print("test accuracy:     %f", testing)
        results.append("decay = " + str(decay) + ", rate = " + str(rate) + ", training accuracy: " + str(training) + ", test accuracy: " + str(testing))

for r in results:
    print(r)

print("Finish")