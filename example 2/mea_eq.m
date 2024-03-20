function z=mea_eq(x)
z=[sqrt(x(1,:).^2+x(2,:).^2);atan2(x(2,:),x(1,:))];
end