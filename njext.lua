D = {}
parent = {}
rogue = {}

function find_closest()
    local mindis = 1231231234
    local minpair = nil
    local R = {}
    local N = 0
    for _, _ in pairs(D) do
        N = N + 1
    end
    local active = {}
    for i, _ in pairs(D) do
        table.insert(active, i)
    end
    for _i = 1,#active do
        for _j = _i+1,#active do
            local i = active[_i]
            local j = active[_j]
            local ip = parent[i]
            local jp = parent[j]
            if rogue[i] or rogue[j] or ip == jp then
                if R[i] == nil then
                    R[i] = 0
                    for k, _ in pairs(D) do
                        R[i] = R[i] + D[i][k]
                    end
                end

                if R[j] == nil then
                    R[j] = 0
                    for k, _ in pairs(D) do
                        R[j] = R[j] + D[j][k]
                    end
                end
                local qij = (N - 2) * D[i][j] - R[i] - R[j]
                if qij < mindis then
                    mindis = qij
                    minpair = {i, j}
                end
            end
        end
    end
    return minpair[1], minpair[2]
end

function join_node(u, v, n)
    for k, _ in pairs(D) do
        D[k][n] = 0.5 * (D[u][k] + D[v][k] - D[u][v])
    end
    D[n] = {}
    for k, _ in pairs(D) do
        if k == n then
            D[n][k] = 0
        end
        D[n][k] = D[k][n]
    end
    -- update parent
    -- if not joined against its direct parent,
    -- then the parent of the new node is n
    if parent[u] ~= n then
        parent[n] = parent[u]
    end
    -- if joined against its direct parent
    -- then no need to update
    D[u] = nil
    D[v] = nil
    local cnt = 0
    for i, _ in pairs(D) do
        cnt = cnt + 1
    end
    return cnt
end