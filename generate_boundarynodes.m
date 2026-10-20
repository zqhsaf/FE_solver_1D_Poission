function boundarynodes = generate_boundarynodes(Pb, left_type, right_type)

n = length(Pb);

boundarynodes = [ left_type,  right_type;
                  1,          n        ];
end
