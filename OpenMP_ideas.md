### OpenMP

#### General Ideas

- First call to `#pragma omp parallel` should be in main when calling `build_tree()`. We don't want to create a new parallel region every time the recursive function is called

```cpp
int main() {
[...]
    #pragma omp parallel
    #pragma omp master  // We don't want all threads to call the function the first time
    build_tree();
[...]
}
```

- Every call to the recursive function inside the body of the function should be a task (`#pragma omp task`)
- `#pragma omp parallel for` to calculate the general performance
  - Discuss if `parallel` is needed or not, might not be a good idea to create parallel sections inside the task, we can run out of threads
- parallelize finding the median and the points allocation for `L` and `R`

#### Experimental Ideas

- Experiment running the provided algorithm to find the 2 furthest points in the set, which is not parallelizable vs running the brute force way, but parallel

```cpp
max = 0.0;
a = null;
b = null;
#pragma omp for 
for(i = 0; i < num_points; i++) {
    for (int j = i; j < num_points; j++) {
        float dist = vec_dist(pts[i], pts[j]);

        #pragma omp critical 
        {
            if (dist > max) {
                max = dist;
                a = pts[i];
                b = pts[j]
            }
        }
    }
}
```
