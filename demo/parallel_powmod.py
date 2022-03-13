import gmpy2
import concurrent.futures
import multiprocessing
import time

# Test parallel computation of gmpy2.powmod()

# Create a set of test values.
def create_tests(num_items = 100000, bits = 1000):
    rand = gmpy2.random_state(42)
    e = gmpy2.mpz_urandomb(rand, bits)
    m = gmpy2.mpz_urandomb(rand, bits)
    threads = multiprocessing.cpu_count()

    big_list = []
    for _ in range(num_items):
        big_list.append(gmpy2.mpz_urandomb(rand, bits))

    return big_list, e, m

# Partition the list into equal-sized sub-lists where each
# sub-list will be allocated to a thread.
def partition_list(big_list, threads):
    remaining_items = len(big_list)
    split_list = []
    offset = 0
    while threads:
        q, r = divmod(remaining_items, threads)
        sublist_size = q + bool(r)
        split_list.append(big_list[offset:offset + sublist_size])
        remaining_items -= sublist_size
        offset += sublist_size
        threads -= 1
    return split_list

# Create the function that will be run in each thread.
def powmod_list(lst, e, m, release_gil = True):
    # Set to True to release the GIL.
    # Set to False to not release the GIL.
    gmpy2.get_context().allow_release_gil = release_gil
    start = time.time()
    for i in lst:
        gmpy2.powmod(i, e, m)
    return time.time() - start

# Create the main function.
# Execute the powmod_list function across multiple threads.
def run_test(big_list, e, m, threads, release_gil = True):
    split_list = partition_list(big_list, threads)
    total_thread_time = 0
    walltime_start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
        tasks = []
        for slist in split_list:
            tasks.append(executor.submit(powmod_list, slist, e, m))
        for task in concurrent.futures.as_completed(tasks):
            total_thread_time += task.result()
    return time.time() - walltime_start, total_thread_time

if __name__ == "__main__":
    big_list, e, m = create_tests(10000, 2000)
    # Specify the number of threads to use.
    # Note: 0 indicates a single thread that does not release the GIL.
    test_threads = [0, 1, 2, 4]
    while test_threads[-1] < 2 * multiprocessing.cpu_count():
        test_threads.append(test_threads[-1] + 4)
    for t in test_threads:
        if t == 0:
            print("Executing reference test without threading or releasing the GIL.")
            print("Elapsed time: ", powmod_list(big_list, e, m, False))
            print()
        else:
            print("Executing tests with threading and releasing the GIL.")
            print("Number of threads: ", t)
            print("Wall time, CPU time: ", run_test(big_list, e, m, t, True))
