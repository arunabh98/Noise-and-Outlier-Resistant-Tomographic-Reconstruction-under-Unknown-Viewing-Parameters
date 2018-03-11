function res = radonTransform(...
    angles, width, height, output_size, projection_length)
    res.angles = angles;
    res.width = width;
    res.height = height;
    res.output_size = output_size;
    res.projection_length = projection_length;
    res.adjoint = 0;

    % Register this variable as a partialDCT class
    res = class(res,'radonTransform');