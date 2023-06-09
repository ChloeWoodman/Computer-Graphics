#pragma once
#include <condition_variable>
#include <functional>
#include <vector>
#include <thread>
#include <queue>

// A class that represents a thread pool
class ThreadPool {
public:
	using Task = std::function<void()>;

	// Constructor to start a thread pool with a given number of threads
	ThreadPool(short threadCount) {
		Start(threadCount);
	}

	// Destructor to stop the thread pool
	~ThreadPool() {
		Stop();
	}

	// Method to add a task to the thread pool queue
	void Enqueue(Task task) {
		{
			std::unique_lock<std::mutex> lock(mMutex);
			mTasks.emplace(std::move(task));
		}

		mCondition.notify_one();
	}

private:
	std::vector<std::thread> mThreads; // The thread pool
	std::condition_variable mCondition; // A condition variable to wait for tasks
	std::mutex mMutex; // A mutex to protect the task queue
	std::queue<Task> mTasks; // A queue to store tasks
	bool mRunning = false; // A flag to stop the threads

    // Method to start the thread pool with a given number of threads
    void Start(short threadCount) {
        for (short i = 0; i < threadCount; i++) {
            // Start a thread and add it to the thread pool
            mThreads.emplace_back([=] {
                while (true) {
                    Task task;
                    {
                        std::unique_lock<std::mutex> lock(mMutex);

                        // Wait for a task or stop signal
                        mCondition.wait(lock, [=] { return mRunning || !mTasks.empty(); });

                        // If there are no tasks and the stop signal is received, break the loop
                        if (mRunning && mTasks.empty())
                            break;

                        // Get the next task from the queue
                        task = std::move(mTasks.front());
                        mTasks.pop();
                    }
                    // Execute the task
                    task();
                }
                });
        }
    }

    // Method to stop the thread pool
    void Stop() {
        {
            std::unique_lock<std::mutex> lock(mMutex);
            mRunning = true;
        }

        // Notify all threads to stop waiting for tasks
        mCondition.notify_all();

        // Wait for all threads to finish
        for (auto& thread : mThreads)
            thread.join();
    }
};
