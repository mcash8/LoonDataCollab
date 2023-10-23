from CommonImports import *

def pruneCapacity(num_nodes, capacity):
    '''
    !!-------For SCF Solver-------!!
    Prune capacity matrix according to how many nodes to investigate 
    Returns a num_nodes x num_nodes capacity matrix 
    '''
    flight_ids = np.load('flight_ids.npy')
    flt_avgs = []
    rows = list(range(0, capacity[0].shape[0]-1)) #don't want to index the last row
    for row in rows: 
        
        temp = capacity[:, row] #get one row through all time indices 
        avg_len = float(len(temp[0])) #length of row for averaging 
        
        #filter row since we're working with the upper triangle of the matrix
        temp = temp[:, (row+1):38] #slicing is non-inclusive so the end idx should be incremented by 1 

        #count number of non-zero elements and average for each row
        temp = np.count_nonzero(temp, axis=1)
        temp = np.divide(temp, avg_len)
        
        #round to 3 decimal places
        temp = np.average(temp)
        temp = np.round(temp, 3)

        #average through time and append to list
        flt_avgs.append(np.average(temp))
    
    
    
    flt_avgs_idx= np.sort(np.argsort(flt_avgs)[-num_nodes:])#index of largest avgs to keep
    flt_ids_reduced = flight_ids[flt_avgs_idx]

    reduced_capacity = capacity[:, flt_avgs_idx, :]
    reduced_capacity = reduced_capacity[:, :, flt_avgs_idx]

    return reduced_capacity, flt_ids_reduced


def getBidirectionalCapacityConstraints(combos, flow_var_dict, topology_var_dict, capacity, solver, infinity):
    '''
    !!-------For SCF Solver-------!!
    Define the topology control and flow variables 
    and add bidirectional and flow constraints to solver
    Returns a dictionary of flow and topology constraints for the solver
    '''

    for combo in combos:
        i, j = combo
        #flow keys
        key0 = 'x' + str(i) +'_'+ str(j) #i,j
        key1 = 'x' + str(j) +'_'+str(i) #j,i

        #add flow variables
        flow_var_dict[key0] = solver.IntVar(0.0, infinity, "x" + str(i)+str(j))
        flow_var_dict[key1] = solver.IntVar(0.0, infinity, "x" + str(j)+str(i))

        #add topology variables 
        key2 = "z" + str(i)+'_'+str(j) #i,j
        key3 = "z" + str(j)+'_'+str(i) #j,i

        topology_var_dict[key2] = solver.IntVar(0.0, 1.0, "z" + str(i)+str(j))
        topology_var_dict[key3] = solver.IntVar(0.0, 1.0, "z" + str(j)+str(i))

        #add bidriectional and capacity constraints 
        solver.Add(topology_var_dict[key2] == topology_var_dict[key3]) #bidirectional constraint

        solver.Add(flow_var_dict[key0] <= capacity[i-1,j-1]*topology_var_dict[key2]) #capacity constraint i,j
        solver.Add(flow_var_dict[key1] <= capacity[i-1,j-1]*topology_var_dict[key3]) #capacity constraint j,i

    #print('Number of variables =', solver.NumVariables())  

    return topology_var_dict, flow_var_dict
    
def getFlowTopologyConstraints(flow_var_dict, topology_var_dict, rho, L_k, solver, T, source, dest, num_nodes):
    '''
    !!-------For SCF Solver-------!!
    Define flow conservation constraint equations and topology control constraints
    Returns a dictionary of flow and topology constraints for the solver
    '''

    #flow conservation constraint equations and topology control constraint
    for i in range(1,num_nodes+1): 
        add = []
        sub = []
        add_top = [] 
        for j in range(1,num_nodes+1): 
            if i!=j: 
                key1 = 'x'+ str(i) + '_' + str(j)
                key2 = 'x' + str(j) + '_' + str(i)

                key3 = 'z' + str(i) +'_'+ str(j)

                add.append(flow_var_dict[key1]) #x_ij
                sub.append(flow_var_dict[key2]) #x_ji

                add_top.append(topology_var_dict[key3]) #z_ij

            else: 
                continue

        solver.Add(sum(add_top) <= T) #topology control

        #flow conservation piecewise function
        if type(source.item())!= int:
            if i in source: 
                solver.Add(sum(add)-sum(sub) == rho*L_k)
            elif i in dest: 
                solver.Add(sum(add)-sum(sub) == -rho*L_k)
            else: 
                solver.Add(sum(add)-sum(sub) == 0)
        else: 
            if i == source: 
                solver.Add(sum(add)-sum(sub) == rho*L_k)
            elif i == dest: 
                solver.Add(sum(add)-sum(sub) == -rho*L_k)
            else: 
                solver.Add(sum(add)-sum(sub) == 0)

    return topology_var_dict, flow_var_dict

def printSolutionValues(num_nodes, flow_var_dict, topology_var_dict, rho):
    '''
    !!-------For Debugging-------!!
    Print solver solution. Enable for debugging purposes. 
    '''
    #print('Solution:')
    #print('p = ', rho.solution_value())
     
    combos = itertools.combinations(range(1,num_nodes+1), 2)
    flow_vals = np.zeros((num_nodes, num_nodes-1))
    topology_vals = np.zeros((num_nodes, num_nodes-1))
   
    for i in range(1,num_nodes+1):
        k=0 #positioning for column, will get messy if we use j variable
        
        for j in range(1, num_nodes+1):
            if i!=j: 
                key1 = 'x' + str(i) + '_' + str(j)
                key2 = 'z' + str(i) + '_' + str(j)
                
                val1 = flow_var_dict[key1] #xij
                val2 = topology_var_dict[key2]
                
                flow_vals[i-1, k] = val1.solution_value()
                topology_vals[i-1, k] = val2.solution_value()
                k+=1
    
    #print('flow values')
    #print(flow_vals)
 
    
    #print('-----------------------------------------------------')

    #print('topology control values')
    #print(topology_vals)
    
    return flow_vals, topology_vals


def addTopologyVariablesAndConstraints(topology_var_dict, T, solver, combos, num_nodes): 
    '''
    !!-------For MCF Solver-------!!
    Define the topology variables and the topology constraints
    Add the topology variables and topology contraints to solver
    '''

    #adding topology variables and bidirectional constraint
    for combo in combos: 
        i, j = combo

        #topology keys 
        key1 = "z" + str(i)+'_'+str(j) #i,j
        key2 = "z" + str(j)+'_'+str(i) #j,i

        topology_var_dict[key1] = solver.IntVar(0.0, 1.0, "z" + str(i)+str(j))
        topology_var_dict[key2] = solver.IntVar(0.0, 1.0, "z" + str(j)+str(i))

        #add bidirectional constraint
        solver.Add(topology_var_dict[key1] == topology_var_dict[key2]) 
    
    #accessing topology variables and adding topology control constraint
    for i in range(1,num_nodes+1): 
        add_top = [] 
        for j in range(1,num_nodes+1): 
            if i!=j: 
                key = 'z' + str(i) +'_'+ str(j)
                add_top.append(topology_var_dict[key]) #z_ij

            else: 
                continue
        #sum over j        
        solver.Add(sum(add_top) <= T) 
    return topology_var_dict

def addFlowVariablesAndConstraints(topology_var_dict, flow_var_dict, rho, L_k, solver, source, dest, combos, infinity, capacity, num_nodes): 
    '''
    !!-------For MCF Solver-------!!
    Define the flow variables and the flow constraints
    Add the flow variables and flow contraints to solver
    '''

    #add flow variables to solver, loop through each commodity, create a nested dictionary for each commodity
    for i in range(L_k):
        outter_key = str(i) 
        inner_dict = {} 
        for combo in combos: 
            j, k = combo
            key1 = 'x'+ str(j) + '_' + str(k)
            key2 = 'x' + str(k) + '_' + str(j)

            #add flow variables
            inner_dict[key1] = solver.IntVar(0.0, infinity, "x" + str(j)+str(k))
            inner_dict[key2] = solver.IntVar(0.0, infinity, "x" + str(k)+str(j))


        #nested dict 
        flow_var_dict[outter_key] = inner_dict
    
    #add capacity constraint
    for combo in combos:
        #for each (i,j) pair
        i, j = combo

        #flow keys
        key0 = 'x' + str(i) +'_'+ str(j) #i,j
        key1 = 'x' + str(j) +'_'+str(i) #j,i

        key2 = "z" + str(i)+'_'+str(j) #i,j
        key3 = "z" + str(j)+'_'+str(i) #j,i

        add_1 = []
        add_2 = [] 

        #sum over k
        for k in range(L_k): 
            inner_key = str(k)
            inner_dict = flow_var_dict[inner_key]

            add_1.append(inner_dict[key0])
            add_2.append(inner_dict[key1])
        
        #capacity constraint
        solver.Add(sum(add_1) <= capacity[i-1,j-1]*topology_var_dict[key2])
        solver.Add(sum(add_2) <= capacity[i-1,j-1]*topology_var_dict[key3])
    

    #add flow conservation constraint
    for k in range(L_k): 
        add = [] 
        sub = [] 

        inner_key = str(k)        
        inner_dict = flow_var_dict[inner_key]

        for i in range(1, num_nodes+1): 

            for j in range(1, num_nodes+1): 
                if i!=j: 
                    key1 = 'x'+ str(i) + '_' + str(j)
                    key2 = 'x' + str(j) + '_' + str(i)

                    add.append(inner_dict[key1]) #x_ij
                    sub.append(inner_dict[key2]) #x_ji
                else: 
                    continue
            
            if i in source: 
                solver.Add(sum(add)-sum(sub) == rho*L_k[k])
            elif i in dest: 
                solver.Add(sum(add)-sum(sub) == -rho*L_k[k])
            else: 
                solver.Add(sum(add)-sum(sub) == 0)

    return topology_var_dict, flow_var_dict

def singleCommodityFlow(num_nodes, source, dest, L_k, T, capacity):
    
    '''
    Solve single commodity flow problem
    '''
    
    solver = pywraplp.Solver.CreateSolver('SAT')
    if not solver:
        return

    infinity = solver.infinity()

    #add rho variable
    rho = solver.NumVar(0.0, infinity, 'p')

    #define flow variables and topology control variables 
    combos = itertools.combinations(range(1,num_nodes+1), 2)
    flow_var_dict = {} 
    topology_var_dict = {} 
 
    #setup the problem
    topology_var_dict, flow_var_dict = getBidirectionalCapacityConstraints(combos, flow_var_dict, topology_var_dict, capacity, solver, infinity)
    
    
    topology_var_dict, flow_var_dict = getFlowTopologyConstraints(flow_var_dict, topology_var_dict, rho, L_k, solver, T, source, dest, num_nodes)
   
         
    #define the objective
    solver.Maximize(rho)

    #call the solver
    #print(f'Solving with {solver.SolverVersion()}')
    status = solver.Solve()
    
    
    if status == pywraplp.Solver.OPTIMAL:
        flow_vals, topology_vals = printSolutionValues(num_nodes, flow_var_dict, topology_var_dict, rho)
 
        '''
        print('\nAdvanced usage:')
        print('Problem solved in %f milliseconds' % solver.wall_time())
        print('Problem solved in %d iterations' % solver.iterations())
        print('Problem solved in %d branch-and-bound nodes' % solver.nodes())
        '''
        return rho.solution_value(), topology_vals
        
    else:
        print('The problem does not have an optimal solution.')
        
        return

def solveMultiCommodityFlow(num_nodes, source, dest, L_k, T, capacity): 
    '''
    Solve multicommodity flow problem 
    '''

    #check source, dest, and l_k are lists 
    if not all(isinstance(i, list) for i in [source, dest, L_k]):
        raise Exception('Source, destination, or commodity variables are not lists')

    #instantiate solver 
    solver = pywraplp.Solver.CreateSolver('SAT')
    if not solver:
        return

    infinity = solver.infinity()

    #add rho variable
    rho = solver.NumVar(0.0, infinity, 'p')

    combos = itertools.combinations(range(1,num_nodes+1), 2) #for creating keys

    #add topology variables and topology constraints 
    topology_var_dict = {} 
    topology_var_dict = addTopologyVariablesAndConstraints(topology_var_dict, T, solver, combos)

    #add flow variables and flow constraints for each commodity
    flow_var_dict = {}

    flow_var_dict = addFlowVariablesAndConstraints(topology_var_dict, flow_var_dict, rho, L_k, solver, source, dest, infinity)

    #define the objective
    solver.Maximize(rho)

    #call the solver
    #print(f'Solving with {solver.SolverVersion()}')
    status = solver.Solve()
    
    
    if status == pywraplp.Solver.OPTIMAL:
        flow_vals, topology_vals = printSolutionValues(num_nodes, flow_var_dict, topology_var_dict, rho)
 
        '''
        print('\nAdvanced usage:')
        print('Problem solved in %f milliseconds' % solver.wall_time())
        print('Problem solved in %d iterations' % solver.iterations())
        print('Problem solved in %d branch-and-bound nodes' % solver.nodes())
        '''
        return rho.solution_value(), topology_vals
        
    else:
        print('The problem does not have an optimal solution.')
        
        return
