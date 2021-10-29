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
v = zeros(n)
tolerance = 0.001
imax = 1000
iter = 1
vnew = zeros(n)
v = vnew .+ 2*tolerance

while maximum(abs.(v - vnew)) > tolerance && iter < imax
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
    end  

    ### STEP 2b
    (vnew, cartesianindex) = findmax(log.(c) .+ β*v', dims = 2) 
    iter += 1 #adding counter
end

vnew

scatter(kgrid, vnew, titke = "V(k)")

#Optimal policy function
kprimeindex = getindex