import gmpy2
import concurrent.futures
import time

# Test parallel computation of gmpy2.powmod()

rand = gmpy2.random_state(42)

num_items = 100_000
bits = 1_000

# Change threads to reflect the number of CPU cores to use.
threads = 16

e = gmpy2.mpz_urandomb(rand, bits)
m = gmpy2.mpz_urandomb(rand, bits)

big_list = []
for _ in range(num_items):
    big_list.append(gmpy2.mpz_urandomb(rand, bits))

# Partition the list into equal-sized sub-lists where each
# sub-list will be allocated to a thread.
split_list = []
temp_threads = threads
temp_num_items = num_items
j = 0
while temp_threads:
    split_size = temp_num_items // temp_threads + bool(temp_num_items % temp_threads)
    split_list.append(big_list[j:j+split_size])
    temp_num_items -= split_size
    j += split_size
    temp_threads -= 1

# Create the function that will be run in each thread.
def powmod_list(lst, e, m):
    # Set to True to release the GIL.
    # Set to False to not release the GIL.
    gmpy2.get_context().allow_release_gil = True
    start = time.time()
    for i in lst:
        gmpy2.powmod(i, e, m)
    return time.time() - start

# Execute the powmod_list function across multiple threads.
with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
    tasks = []
    for slist in split_list:
        tasks.append(executor.submit(powmod_list, slist, e, m))
    for task in concurrent.futures.as_completed(tasks):
        print(task.result())
