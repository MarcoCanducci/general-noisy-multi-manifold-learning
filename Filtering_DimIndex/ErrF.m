function F = ErrF

v1 = [1 0 0];
v2 = [0.5 0.5 0];
v3 = [1/3 1/3 1/3];

F = @(x)(pdist2(x,v1,@FisherMetricDist) - pdist2(x,v2,@FisherMetricDist))^2 + ...
        (pdist2(x,v2,@FisherMetricDist) - pdist2(x,v3,@FisherMetricDist))^2 + ...
        (pdist2(x,v1,@FisherMetricDist) - pdist2(x,v3,@FisherMetricDist))^2;

