""" Neural Network.
A 2-Hidden Layers Fully Connected Neural Network (a.k.a Multilayer Perceptron)
implementation with TensorFlow. This example is using the MNIST database
of handwritten digits (http://yann.lecun.com/exdb/mnist/).
Links:
    [MNIST Dataset](http://yann.lecun.com/exdb/mnist/).
Author: Aymeric Damien
Project: https://github.com/aymericdamien/TensorFlow-Examples/
"""
def get_batch(total, batch_number):
    sample = np.arange(total)
    np.random.shuffle(sample)
    batch = np.array_split(sample,batch_number)
    return batch


from __future__ import print_function

# Import MNIST data
from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets("/tmp/data/", one_hot=True)

import tensorflow as tf
import numpy as np

# Parameters
learning_rate = 0.1
num_steps = 1000
batch_size = 128
display_step = 100

# Network Parameters
n_hidden_1 = 256 # 1st layer number of neurons
n_hidden_2 = 256 # 2nd layer number of neurons
num_input = 784 # MNIST data input (img shape: 28*28)
num_classes = 10 # MNIST total classes (0-9 digits)

# tf Graph input
X = tf.placeholder("float", [None, num_input])
Y = tf.placeholder("float", [None, num_classes])

# Store layers weight & bias
# weights = {
#     'h1': tf.Variable(tf.random_normal([num_input, n_hidden_1])),
#     'h2': tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2])),
#     'out': tf.Variable(tf.random_normal([n_hidden_2, num_classes]))
# }
# biases = {
#     'b1': tf.Variable(tf.random_normal([n_hidden_1])),
#     'b2': tf.Variable(tf.random_normal([n_hidden_2])),
#     'out': tf.Variable(tf.random_normal([num_classes]))
# }


# Create model
def neural_net(x,lam,w2):
    # Hidden fully connected layer with 256 neurons
    layer_1 = tf.layers.dense(x,128,activation=tf.nn.relu,
                              kernel_regularizer=tf.contrib.layers.l1_regularizer(scale=lam * w2))
    # Hidden fully connected layer with 256 neurons
    layer_2 = tf.layers.dense(layer_1,64,activation=tf.nn.relu,
                              kernel_regularizer=tf.contrib.layers.l1_regularizer(scale=lam * w2))

    # Output fully connected layer with a neuron for each class
    out_layer = tf.layers.dense(layer_2,num_classes,activation=None,
                                kernel_regularizer=tf.contrib.layers.l1_regularizer(scale=lam*w2))
    return out_layer

# Construct model


train_name = "train_1.txt"
test_name = "test_1.txt"
train_acc = []
test_acc = []
with tf.Session() as sess:
    for t in range(100):
        print(t)
        weight = np.random.exponential(size=(60000, 1))
        weight1 = weight / np.mean(weight)
        weight2 = np.random.exponential(size=1)
        logits = neural_net(X, lam=0.0001, w2=weight2)
        prediction = tf.nn.softmax(logits)

        # Define loss and optimizer
        W1 = tf.placeholder("float", [None, 1])
        loss_op = tf.reduce_mean(tf.multiply(tf.nn.softmax_cross_entropy_with_logits(
            logits=logits, labels=Y), W1))
        optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate)
        train_op = optimizer.minimize(loss_op)

        # Evaluate model
        correct_pred = tf.equal(tf.argmax(prediction, 1), tf.argmax(Y, 1))
        accuracy = tf.reduce_mean(tf.cast(correct_pred, tf.float32))

        # Initialize the variables (i.e. assign their default value)
        init = tf.global_variables_initializer()

        # get batch
        batch_number = int(60000 / batch_size)
        batch = get_batch(60000, batch_number)

        # Start training

        # Run the initializer
        sess.run(init)
        idx = 0

        for step in range(1, num_steps + 1):
            batch_x, batch_y = mnist.train.next_batch(batch_size)
            # Run optimization op (backprop)
            sess.run(train_op, feed_dict={X: batch_x, Y: batch_y, W1: weight1[batch[idx]]})
            idx += 1
            idx = idx % batch_number
            if step % display_step == 0 or step == 1:
                # Calculate batch loss and accuracy
                loss, acc = sess.run([loss_op, accuracy], feed_dict={X: batch_x,
                                                                     Y: batch_y,
                                                                     W1: weight1[batch[idx]]})
                print("Step " + str(step) + ", Minibatch Loss= " + \
                      "{:.4f}".format(loss) + ", Training Accuracy= " + \
                      "{:.3f}".format(acc))

        print("Optimization Finished!")
        train = sess.run(accuracy, feed_dict={X: mnist.train.images, Y: mnist.train.labels})
        test = sess.run(accuracy, feed_dict={X: mnist.test.images, Y: mnist.test.labels})
        train_acc.append(train)
        test_acc.append(test)

        # Calculate accuracy for MNIST test images
        print("Testing Accuracy:", test)

np.savetxt(train_name, train_acc)
np.savetxt(test_name, test_acc)




