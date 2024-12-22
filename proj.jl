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
    if isnothing(death_sx)
        println("Bar is infinite")
        return
    end

    simplices = cloud[death_sx]

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
       
        cloud[findall(x -> x == simplices[i], cloud)[1]] = (simplices[i][1] - avg[1]*epsilon*direction, simplices[i][2] - avg[2]*epsilon*direction)
        #println(cloud[death_sx][i])
        #println("................")
    end
    println(cloud)
end


function influence_birth(bar, cloud, epsilon = 0.1, direction = 1)
    #direction -1 shortens the bar by making it die sooner, 1 lengthens the bar by making it die later

    birth_simplex = PersistenceDiagrams.birth_simplex(bar)
    println(birth_simplex)
    if length(birth_simplex) == 1
        #the birth simplex is a single pint in 1D homology, they all start at 0
        return
    end
    
    print(cloud[birth_simplex])

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
        println(avg)
        cloud[findall(x -> x == simplices[i], cloud)[1]] = (simplices[i][1] - avg[1]*epsilon*direction, simplices[i][2] - avg[2]*epsilon*direction)

    end

    
end




cloud = random_cloud(10)
#println(cloud)
scatter(cloud)
pd = bars(cloud)
most_persistent = pd[2][end]
#println(most_persistent)
#println("test0")
influence_death(most_persistent, cloud)

#println("TEST")
#println(cloud)
scatter!(cloud)

influence_birth(most_persistent, cloud)

scatter!(cloud)


