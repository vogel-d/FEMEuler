function insertFaceFace!(incfe::Array{Int,1},EdgeNumberStart1::Int,EdgeNumberStart2::Int, EdgeNumberStartEW1::Int, EdgeNumberStartEW2::Int,
                         EdgeNumberStartSN1::Int, EdgeNumberStartSN2::Int, n::Int)

  EdgeNumber1=EdgeNumberStart1
  EdgeNumber2=EdgeNumberStart2
  EdgeNumberEW1=EdgeNumberStartEW1
  EdgeNumberEW2=EdgeNumberStartEW2
  EdgeNumberSN1=EdgeNumberStartSN1
  EdgeNumberSN2=EdgeNumberStartSN2
  for j in 1:n
    for i in 1:n
      #=
      if ((i==1 && j==1) || (i==1 && j==n) || (i==n && j==1) || (i==n && j==n))
        Face(FaceNumber)%CornerCube=.TRUE.
      end
      =#
      #Face(FaceNumber)%Number=FaceNumber
      #Face(FaceNumber)%Boundary=1

      e=Int[EdgeNumber1-n, EdgeNumber1, EdgeNumber2-1, EdgeNumber2];
      if j==1
        e[1]=EdgeNumberEW1
        EdgeNumberEW1+=1
      end
      if j==n
        e[2]=EdgeNumberEW2
        EdgeNumberEW2+=1
      end
      if i==1
        e[3]=EdgeNumberSN1
        EdgeNumberSN1+=1
      end
      if i==n
        e[4]=EdgeNumberSN2
        EdgeNumberSN2+=1
        EdgeNumber2-=1
      end
      append!(incfe,e)
      EdgeNumber1+=1
      EdgeNumber2+=1
    end
  end
  return nothing;
end

function insertFaceEdge!(ince::Array{Int,1},NodeNumberStart::Int, NodeNumberE1Start1::Int, NodeNumberE2Start1::Int, NodeNumberE1Start2::Int,
                        NodeNumberE2Start2::Int, EdgeNumber::Int, n::Int)

  NodeNumber=NodeNumberStart
  NodeNumberE1=NodeNumberE1Start1
  NodeNumberE2=NodeNumberE2Start1
  EdgeNumberStart1=EdgeNumber;
  for j in 1:n-1
    for i in 1:n
      e=Int[NodeNumber-1, NodeNumber]
      if i==1
        e[1]=NodeNumberE1
        NodeNumberE1+=1
      end
      if i==n
        e[2]=NodeNumberE2
        NodeNumberE2+=1
        NodeNumber-=1
      end
      append!(ince,e)
      EdgeNumber+=1
      NodeNumber+=1
    end
  end
  NodeNumber=NodeNumberStart
  NodeNumberE1=NodeNumberE1Start2
  NodeNumberE2=NodeNumberE2Start2
  EdgeNumberStart2=EdgeNumber;
  for j in 1:n
    for i in 1:n-1
      e=Int[NodeNumber-(n-1), NodeNumber]
      if j==1
        e[1]=NodeNumberE1
        NodeNumberE1+=1
      end
      if j==n
        e[2]=NodeNumberE2
        NodeNumberE2+=1
      end
      append!(ince,e)
      EdgeNumber+=1
      NodeNumber+=1
    end
  end
  return EdgeNumber, EdgeNumberStart1, EdgeNumberStart2
end

function insertEdgeEdge!(ince::Array{Int,1},NodeNumberStart::Int, NodeNumberE1::Int, NodeNumberE2::Int, EdgeNumber::Int, n::Int)
  NodeNumber=NodeNumberStart
  EdgeNumberStart=EdgeNumber
  for i in 1:n
    e=Int[NodeNumber-1, NodeNumber]
    if i==1
      e[1]=NodeNumberE1
    end
    if i==n
      e[2]=NodeNumberE2
    end
    append!(ince,e)
    EdgeNumber+=1
    NodeNumber+=1
  end
  return EdgeNumber, EdgeNumberStart
end
