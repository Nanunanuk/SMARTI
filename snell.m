function vector_out=snell(layer_data,layer_position,adj_layer_position,wavelength,vector_in,surface_vector)
% If the adjajent medium is the ambient type 0 for adj_layer_position

if adj_layer_position>length(layer_data)
    n1=layer_data(layer_position).layer_mat(find(layer_data(layer_position).layer_mat==wavelength,1),2);
    n2=1;
    
elseif adj_layer_position==0    
    n1=layer_data(layer_position).layer_mat(find(layer_data(layer_position).layer_mat==wavelength,1),2);
    n2=1;
    
else
    n1=layer_data(layer_position).layer_mat(find(layer_data(layer_position).layer_mat==wavelength,1),2);
    n2=layer_data(adj_layer_position).layer_mat(find(layer_data(adj_layer_position).layer_mat==wavelength,1),2);
end

d=vector_in;
d=d/norm(d);
n=surface_vector;
n=n/norm(n);

vector_out=n1/n2 *(d+n*dot(-n,d))-n*sqrt(1-(n1/n2)^2*(1-dot(-n,d)^2));      % (Marschner,Shirley,Fundamentals of Computer Graphics,2016) 


