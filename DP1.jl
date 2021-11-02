###########################
### DYNAMIC PROGRAMMING ###
###########################


### STEP 0
α = 0.3
δ = 0.1
β = 0.9

kupper = 2
klower = 0.001
n = 10
kgrid = collect(range(klower, stop=kupper, length=n))

### STEP 1 
tolerance = 0.001
imax = 1000
vnew = zeros(length(kgrid)) # initial value function guess
v = vnew .+ 2*tolerance # initialize v in a way that ensures loop will start
global cartesianindex = Array{CartesianIndex{2}, n}
i = 1

while maximum(abs.(v - vnew)) > tolerance && i <= imax
   ### STEP 2a - calculate frist by columns in the first row and then jump to the other rwo
    # pre-allocate memory
    v = vnew
    c = zeros(n, n) 
    for i in 1:n
        for j in 1:n
            c[i,j] = kgrid[i]^α + (1-δ)*kgrid[i] - kgrid[j]
            if c[i,j] < 0
                c[i,j] = 0
            end    
        end
    ### STEP 2b
    (vnew, cartesianindex) = findmax(log.(c) .+ β*v', dims = 2) 
    end  
    i += 1 #adding counter
end

vnew

scatter(kgrid, vnew, title = "v(k)")

# Policy function
kprimeindex = getindex.(cartesianindex, 2)
kprime = kgrid[kprimeindex]
scatter(kgrid, kprime, title = "k'(k)")