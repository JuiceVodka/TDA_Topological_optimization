using LinearAlgebra, Plots, Ripserer, PersistenceDiagrams

function random_cloud(n)
    cloud = Tuple{Float64, Float64}[]
    for i in 1:n
        push!(cloud, Tuple(randn(2)))
    end

    return cloud
end

function bars(cloud)
    pd = ripserer(cloud)
    #for i in 1:length(pd[1])
    #    println(pd[1][i])
    #end
    return pd
end


function influence_death(bar, cloud, epsilon = 0.1, direction = 1)
    #direction -1 shortens the bar by making it die sooner, 1 lengthens the bar by making it die later

    #to lengthen a bar, we need to move the faces of its death simplex further apart
    death_sx = death_simplex(bar)
    birth_simplex = PersistenceDiagrams.birth_simplex(bar)
    #print(birth_simplex, " ", death_sx)
    if isnothing(death_sx)
        println("Bar is infinite")
        return
    end

    simplices = cloud[death_sx]
    birth_simplices = cloud[birth_simplex]
    #println(simplices, " ", birth_simplices)

    #move the simplices apart from eachother
    for i in 1:length(simplices)
        vectors = Tuple{Float64, Float64}[]
        for j in 1:length(simplices)
            if i != j
                push!(vectors, (simplices[j][1] - simplices[i][1], simplices[j][2] - simplices[i][2]))
            end
        end
        #find the average vector
        avg = (0, 0)
        for j in 1:length(vectors)
            avg = (avg[1] + vectors[j][1], avg[2] + vectors[j][2])
        end
        avg = (avg[1]/length(vectors), avg[2]/length(vectors))
        #move the point in the direction of the average vector
        #println(avg)
        #println(cloud[death_sx][i])
        #println(simplices[i])
        #println((simplices[i][1] + avg[1]*epsilon, simplices[i][2] + avg[2]*epsilon))
        #println(simplices[i], " ", birth_simplices[1], " ", birth_simplices[2])
        #if length(birth_simplices) == 2 && simplices[i] ∈ birth_simplices
        #    #println("FOUND")
        #    continue
        #end
        #println(simplices[i], " ", birth_simplices[1], " ", birth_simplices[2])

        cloud[findall(x -> x == simplices[i], cloud)[1]] = (simplices[i][1] - avg[1]*epsilon*direction, simplices[i][2] - avg[2]*epsilon*direction)
        #println(cloud[death_sx][i])
        #println("................")
    end
    # println(cloud)
end


function influence_birth(bar, cloud, epsilon = 0.1, direction = 1)
    #direction -1 shortens the bar by making it die sooner, 1 lengthens the bar by making it die later

    birth_simplex = PersistenceDiagrams.birth_simplex(bar)
    # println(birth_simplex)
    if length(birth_simplex) == 1
        #the birth simplex is a single pint in 1D homology, they all start at 0
        return
    end
    
    # print(cloud[birth_simplex])

    simplices = cloud[birth_simplex]

    #move the simplices apart from eachother to make it born later or closer to make it born sooner
    for i in 1:length(simplices)
        vectors = Tuple{Float64, Float64}[]
        for j in 1:length(simplices)
            if i != j
                push!(vectors, (simplices[j][1] - simplices[i][1], simplices[j][2] - simplices[i][2]))
            end
        end
        #find the average vector
        avg = (0, 0)
        for j in 1:length(vectors)
            avg = (avg[1] + vectors[j][1], avg[2] + vectors[j][2])
        end
        avg = (avg[1]/length(vectors), avg[2]/length(vectors))
        # println(avg)
        cloud[findall(x -> x == simplices[i], cloud)[1]] = (simplices[i][1] - avg[1]*epsilon*direction, simplices[i][2] - avg[2]*epsilon*direction)

    end

    
end


function circular_pointcloud(n, r)
    cloud  = Tuple{Float64, Float64}[]
    for i in 1:n
        theta = 2*pi*rand()
        push!(cloud, (r*cos(theta), r*sin(theta)))
    end
    return cloud
end

function circle_bounded_pointcloud(n, r)
    cloud  = Tuple{Float64, Float64}[]
    for i in 1:n
        theta = 2*pi*rand()
        R = r*rand()
        push!(cloud, (R*cos(theta), R*sin(theta)))
        
    end
    return cloud
end


# means birth and death are the same
function is_diagonal(q)
    return abs(q[1] - q[2]) < 1e-6  # Small tolerance for floating-point comparison
end


function adjust_points(A, B) 
    a0, a1 = ripserer(A)
    iterations = 1000
    tolerance = 0.01
    epsilon = 0.1
    for i in 1:iterations
        b0, b1 = ripserer(B)
        
        println("iteration $(i)")
        # Compute matchings for both 0-dim and 1-dim barcodes
        match_bottle_0 = matching(Bottleneck(), a0, b0)
        match_bottle_1 = matching(Bottleneck(), a1, b1)

        for (p, q) in matching(match_bottle_0)
            
            if (is_diagonal(p)) && (is_diagonal(q))
                continue
            end


            # Compare births
            if q.birth < p.birth
                influence_birth(q, B, epsilon, 1)  # Delay birth
            elseif q.birth > p.birth
                influence_birth(q, B, epsilon, -1)  # Advance birth
            end

            # Compare deaths
            if q.death < p.death
                influence_death(q, B, epsilon, 1)  # Prolong death
            elseif q.death > p.death
                influence_death(q, B, epsilon, -1)  # Shorten death
            end
       
        end

        for (p, q) in matching(match_bottle_1)

            if (is_diagonal(p)) && (is_diagonal(q))
                continue
            end

            # Compare births
            if q.birth < p.birth
                influence_birth(q, B, epsilon, 1)  # Delay birth
            elseif q.birth > p.birth
                influence_birth(q, B, epsilon, -1)  # Advance birth
            end

            # Compare deaths
            if q.death < p.death
                influence_death(q, B, epsilon, 1)  # Prolong death
            elseif q.death > p.death
                influence_death(q, B, epsilon, -1)  # Shorten death
            end
        end

        # # Compute new distance
        # match_bottle_0 = matching(Bottleneck(), a0, b0)
        loss = weight(match_bottle_0)
        println("Current loss: $loss")
    
        # Stop if loss is below tolerance
        if loss < tolerance
            break
        end
    end
end    



r = 1
n = 100
A = circular_pointcloud(n, r)
B = circle_bounded_pointcloud(n, r)

scatter(A)
scatter!(B)

#adjust_points(A, B)



# iterations = 100
# tolerance = 0.01
# for _ in 1:iterations
#     b0, b1 = ripserer(B)  

#     match_bottle = matching(Bottleneck(), a0, b0)
#     loss = distance(match_bottle)
#     println("Current loss: $loss")

#     # Stop if the loss is below the tolerance
#     if loss < tolerance
#         break
#     end
# end

#scatter!(B, label="B")

# Display PDs
# pd = bars(A)
# plot(pd, label="A")
# pd = bars(B)
# plot!(pd, label="B")







# scatter(A)
# scatter!(B)


cloud = random_cloud(10)
# #println(cloud)
# scatter(cloud)
pd = bars(cloud)
most_persistent = pd[2][end]
# #println(most_persistent)
# #println("test0")
influence_death(most_persistent, cloud)

# #println("TEST")
# #println(cloud)
# scatter!(cloud)

# influence_birth(most_persistent, cloud)

# scatter!(cloud)