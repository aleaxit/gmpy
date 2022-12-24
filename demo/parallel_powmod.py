import gmpy2
import concurrent.futures
import time

# Test parallel computation of gmpy2.powmod_base_list()

# The demo code is not (yet) optimal since it doesn't take advantage of
# read-only access to data list that is shared across all threads.

# Create a set of test values.
def create_tests(num_items = 100000, bits = 1000):
    rand = gmpy2.random_state(42)
    e = gmpy2.mpz_urandomb(rand, bits)
    m = gmpy2.mpz_urandomb(rand, bits)
    big_list = [gmpy2.mpz_urandomb(rand, bits) for _ in range(num_items)]
    return big_list, e, m

# Partition big_list into a new list containing a number of sub-list with
# maximum length of size.
def partition_list(big_list, size):
    remaining_items = len(big_list)
    split_list = []
    offset = 0
    while offset < len(big_list):
        split_list.append(big_list[offset:offset + size])
        offset += size
    return split_list

# This is the non-threaded, baseline function that does not
# the GIL. It is intended to use the original, non-partitioned list.
def powmod_list_gil(lst, e, m):
    gmpy2.get_context().allow_release_gil = False
    result = []
    start = time.time()
    result = [gmpy2.powmod(i, e, m) for i in lst]
    return time.time() - start, result

# This is the threaded, baseline function that does release the GIL.
def powmod_list_nogil(index, lst, e, m):
    gmpy2.get_context().allow_release_gil = True
    result = []
    start = time.time()
    result = [gmpy2.powmod(i, e, m) for i in lst]
    return time.time() - start, index, result

# This function uses the vector version powmod_base_list and releases the GIL.
def powmod_vector_nogil(index, vector, e, m):
    gmpy2.get_context().allow_release_gil = True
    start = time.time()
    result = gmpy2.powmod_base_list(vector, e, m)
    return time.time() - start, index, result

# Run threaded versions.
def run_test(function, big_list, e, m, threads, release_gil = True):
    over_split = 8
    size = len(big_list) // (over_split * threads) + bool(len(big_list) % (over_split * threads))
    split_list = partition_list(big_list, size)
    total_thread_time = 0
    walltime_start = time.time()
    executor = concurrent.futures.ThreadPoolExecutor(max_workers = threads)
    with executor:
        tasks = []
        for index, slist in enumerate(split_list):
            tasks.append(executor.submit(function, index, slist, e, m))
        for task in concurrent.futures.as_completed(tasks):
            elapsed, index, result = task.result()
            total_thread_time += elapsed
    return time.time() - walltime_start, total_thread_time

if __name__ == "__main__":
    big_list, e, m = create_tests(10000, 1024)

    # Run a baseline test with releasing the GIL.
    # print("Baseline test without releasing the GIL. ", powmod_list_gil(big_list, e, m)[0])

    # Specify the number of threads to use.
    test_threads = [1,2,4,8,12,16]
    # The following demo code is disabled by default. It high-lights the inefficiency of 
    # releasing the GIL for functions that execute fairly quicky
    for i in range(0):
        print("Threaded, list-based, releasing the GIL, pass: ", i + 1)
        for t in test_threads:
            # print("Executing tests with threading and releasing the GIL.")
            print("Number of threads: ", t, end = '')
            print(" Wall time, CPU time: ", run_test(powmod_list_nogil, big_list, e, m, t, True))

    # Repeat the tests multiple times to try to trigger a crash.
    for i in range(100):
        print("Threaded, vector-based, releasing the GIL, pass: ", i + 1)
        for t in test_threads:
            # print("Executing tests with threading and releasing the GIL.")
            print("Number of threads: ", t, end = '')
            print(" Wall time, CPU time: ", run_test(powmod_vector_nogil, big_list, e, m, t, True))
