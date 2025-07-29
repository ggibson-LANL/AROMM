def triangular_lhs(num_samples, bounds):
    """
    Generate random samples from a triangular distribution using Latin Hypercube Sampling (LHS).

    Parameters:
    - num_samples: Number of samples to generate.
    - bounds: A list of tuples, where each tuple contains the lower and upper bounds for each parameter.

    Returns:
    - An array representing the Latin Hypercube Samples from a triangular distribution.
    """
    num_parameters = len(bounds)
    samples = lhs(num_parameters, samples=num_samples)
    #print(samples)
    #print(samples.shape)
  
    
    # Transform LHS samples to triangular distribution
    sum = []
    for i in range(num_parameters):
        samples[:, i] = triang.ppf(samples[:, i], c=(bounds[i][2] - bounds[i][0]) / (bounds[i][1] - bounds[i][0]),
                                   loc=bounds[i][0], scale=bounds[i][1] - bounds[i][0])
        print(samples[:, i],i)
        
        #print(samples[:, i].shape)

        #print(samples[0, i],i)
        
    return samples


def adjust_samples_to_sum(samples, target_sum=1):
       
    while not np.all(success):
        for i in range(num_samples):
            if success[i]:
                
                continue

            def constraint(x):
                print("Start",np.sum(x),"End")
                return np.sum(x) - target_sum

            cons = {'type': 'eq', 'fun': constraint}
            bounds = [(0, 1) for _ in range(num_parameters)]

            result = minimize(
                lambda x: np.sum((samples[i] - x) ** 2),
                samples[i],
                method='SLSQP',
                bounds=bounds,
                constraints=cons
            )
            adjusted_samples[i] = result.x
            if np.isclose(np.sum(result.x), target_sum):
                success[i] = True

        if not np.all(success):
            samples = triangular_lhs(num_samples, bounds)

    return adjusted_samples