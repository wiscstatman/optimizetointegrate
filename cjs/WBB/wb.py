import numpy as np
import tensorflow as tf


def get_batch(total, batch_number):
    sample = np.arange(total)
    np.random.shuffle(sample)
    batch = np.array_split(sample,batch_number)
    return batch


def wbb(xdata, ydata, l, t, epoch,rate):

    n,p = xdata.shape

    x = tf.placeholder(tf.float32, [None, p])
    y = tf.placeholder(tf.float32, [None, 1])
    w = tf.placeholder(tf.float32,[None, 1])
    wp = tf.placeholder(tf.float32)
    beta = tf.Variable(tf.random_normal([p,1], stddev=0.1))
    yhat = tf.matmul(x,beta)
    loss1 = tf.losses.mean_squared_error(y,yhat)
    loss2 = l*wp*tf.reduce_sum(tf.square(beta))
    loss = loss1 + loss2

    batch_number = int(n / 50)
    train = tf.train.GradientDescentOptimizer(rate).minimize(loss)

    beta_estimation = np.zeros((t, p))

    with tf.Session() as sess:
        for i in range(t):
            print(t)
            sess.run(tf.global_variables_initializer())
            weight = np.random.exponential(size=(n,1))
            weight1 = weight / np.mean(weight)
            weight2 = np.random.exponential(size=1)

            for j in range(epoch):
                batch = get_batch(n, batch_number)
                for idx in range(batch_number):
                    _, likelihood, prior = sess.run([train, loss1, loss2],
                                                    feed_dict={x: xdata[batch[idx]], y: ydata[batch[idx]],
                                                               w: weight1[batch[idx]],wp: weight2})
                    print(likelihood + prior)

            current_beta = sess.run(beta)
            beta_estimation[i] =  np.squeeze(current_beta)

    return beta_estimation


data = np.loadtxt("newdata.txt")
xdata = data[:,0:10]
ydata = np.expand_dims(data[:,10],1)


beta = wbb(xdata,ydata,l=1,t=1,epoch=20,rate=0.5)





