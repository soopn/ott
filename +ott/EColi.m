classdef EColi
    % TODO: custom colourmap that has the 0 as the median color?
    properties
        T % T-matrix representation of the particle
    end

    methods
        function this = EColi(wavelength0, n_medium, argv)
        %ECOLI - constructs an E.Coli particle instance

            arguments 
                wavelength0 (1,1) double
                n_medium (1,1) double
                
                argv.n_particle (1,1) double = 1.388
                argv.n_pts (1,1) single = 12
                argv.length (1,1) double = 2
                argv.radius (1,1) double = 1
                argv.Nmax (1,1) single = 14
                argv.method {mustBeText} = 'ebcm'
            end

            scale = 1.0e-6;

            [rho, z] = EColi.parameterise_shape(argv.n_pts, argv.length, argv.radius);

            params = [ rho(:) z(:) ] * scale;
            shape = ott.shapes.AxisymLerp(params);

            this.T = ott.Tmatrix.simple(shape, ...
                'wavelength0', wavelength0, ...
                'index_medium', n_medium, ...
                'index_particle', argv.n_particle, ...
                'method', argv.method, ...
                'Nmax', argv.Nmax ...
            );
        end
    end

    methods (Static) 
        function T = T_matrix(wavelength0, n_medium, argv)
        %T_MATRIX - generates a T-matrix representation of an E.Coli particle
        % Input Arguments:
        % - wavelenght0 (double) -
        %   Initial wavelength of the input beam
        %
        % - n_medium (double) -
        %   Refractive index of the medium
        %
        % - options (name-value-pairs) - 
        %   Available options:
        %   - n_particle (double) - 
        %     Refractive index of the E.Coli particle
        %
        %   - n_pts (single) - 
        %     Number of points to parameterise particle, recommended that
        %     n_pts be an even number greater than 8
        %
        %   - length (double) -
        %     Length of the cylindrical body of the particle in 
        %     unitless dimensions
        %   
        %   - radius (double) -
        %     Radius of cylindrical body and end caps of the particle in
        %     unitlesss dimensions
        %
        %   - Nmax (single) -
        %     Nmax to truncate T-matrix 
        %
        %   - method (string) - 
        %     Method to initialise T-matrix, works best with 'ebcm' or 'pm'
        %
        % Output Arguments:
        % - T (T_matrix) - A Tmatrix object specified by the 'method' input
        
            arguments 
                wavelength0 (1,1) double
                n_medium (1,1) double
                
                % optional parameters
                argv.n_particle (1,1) double = 1.388
                argv.n_pts (1,1) single = 12
                argv.length (1,1) double = 2
                argv.radius (1,1) double = 1
                argv.Nmax (1,1) single = 14
                argv.method {mustBeText} = 'ebcm'
            end

            scale = 1.0e-6;

            [rho, z] = EColi.parameterise_shape(argv.n_pts, argv.radius, argv.length);

            params = [ rho(:) z(:) ] * scale;
            shape = ott.shapes.AxisymLerp(params);

            T = ott.Tmatrix.simple(shape, ...
                'wavelength0', wavelength0, ...
                'index_medium', n_medium, ...
                'index_particle', argv.n_particle, ...
                'method', argv.method, ...
                'Nmax', argv.Nmax ...
            );
        end

        function plots = fixed_beam(beam, wavelength0, n_medium, n_figs, argv)
        %FIXED_BEAM - Generates plots using a fixed input beam 
        %
        % Input Arguments:
        % - beam (Bsc) -
        %   Input beam 
        %
        % - wavelength0 (double) -
        %   Initial wavelength0 of the input beam 
        % 
        % - n_medium (double) -
        %   Refractive index of the medium
        %
        % - n_figs (single) -
        %   Number of desired figures to be generated 
        %
        % - options (name-value pairs) - 
        %   Available options:
        %   
        %   - log_scale (logical) -
        %     Produces output images with a log scale applied
        %
        %   - min_len (double) - 
        %     Minimum length of the E.Coli particle 
        %
        %   - max_len (double) - 
        %     Maximum length of the E.Coli particle
        %
        %   - mean_len (double) - 
        %     Average length of the E.Coli particle
        %
        %   - stdev_len (double) -
        %     Standard deviation for the length of the E.Coli particle
        %
        %   - min_rad (double) - 
        %     Minimum radius of the E.Coli particle 
        %
        %   - max_rad (double) - 
        %     Maximum radius of the E.Coli particle
        %
        %   - mean_rad (double) - 
        %     Average radius of the E.Coli particle
        %
        %   - stdev_rad (double) -
        %     Standard deviation for the radius of the E.Coli particle
        %
        %   - beam_distance (double) -
        %     Distance of the beam from the particle 
        %
        %   - axis (string) -
        %     Axis to view the beam from, specified only by 'x', 'y', or, 'z'
        %
        %   - range (1,2 double) -
        %     Range of values to zoom the image in on
        %
        %   - seed (double) -
        %     Seed to have deterministic randomising of values

            arguments 
                beam ott.Bsc
                wavelength0 double 
                n_medium double
                n_figs single 

                % optional parameters
                argv.log_scale logical = false
                argv.min_len double = 1.0
                argv.max_len double = 2.0
                argv.mean_len double = 4.0
                argv.stdev_len double = 0.08
                argv.min_rad double = 0.35
                argv.max_rad double = 1.0
                argv.mean_rad double = 0.89
                argv.stdev_rad double = 0.08
                argv.beam_distance (1,1) double = 300
                argv.axis {mustBeText} = 'y'
                argv.range (1,2) double = NaN()
                argv.seed double = 0
            end

            rng(argv.seed, "twister")

            % normal distirbution of len and rad
            len = argv.mean_len + argv.stdev_len .* randn(1, n_figs);
            rad = argv.mean_rad + argv.stdev_rad .* randn(1, n_figs);

            % allocate plots
            plots = gobjects(1, n_figs);
            masterFig = figure("Visible", "on");
            master = tiledlayout(masterFig, "flow", ...
                "TileSpacing","compact", ...
                "Padding","compact");
            colormap("gray")

            for i = 1:n_figs

                T = EColi.T_matrix(wavelength0, n_medium, ...
                    "length", len(i), "radius", rad(i));

                sbeam = T * beam;

                % Array of singular figures
                fig = figure('Visible','off');
                colormap("gray")
                plots(i) = fig;

                % Plot fields
                EColi.plot_scatter(fig, sbeam, argv.axis, argv.range, argv.log_scale);
                EColi.plot_total(fig, sbeam, argv.axis, argv.range, argv.log_scale)

                plotParams = ...
                    sprintf('Parameters: $\\ell = %0.2f,~ r = %0.2f$', len(i), rad(i));
                sgtitle(fig, plotParams, 'Interpreter', 'latex')
                

                % Tiled Figures
                tileFig = nexttile(master);

                sbeam.basis = 'regular';
                if isnan(argv.range)
                    sbeam.visualise('axis', argv.axis, 'log_scale', argv.log_scale);
                else
                    sbeam.visualise('axis', argv.axis, 'range', argv.range, 'log_scale', argv.log_scale);
                end

                caption = sprintf(...
                    '$\\ell = %0.2f,~ r = %0.2f$', ...
                    len(i), rad(i) ...
                );
                title(tileFig, caption, 'Interpreter', 'latex'); 
                colorbar
            end
        end

        function plots = image_particle(beam, tMatrix, n_figs, argv)
        %IMAGE_PARTICLE - Generates plots using a given T-matrix
        %
        % Input Arguments:
        % - beam (Bsc) -
        %   Input beam 
        %
        % - wavelength0 (double) -
        %   Initial wavelength0 of the input beam 
        % 
        % - tMatrix (T_matrix) -
        %   Input T-matrix to image
        %
        % - n_figs (single) -
        %   Number of desired figures to be generated 
        %
        % - options (name-value pairs) - 
        %   Available options:
        %   
        %   - log_scale (logical) -
        %     Produces output images with a log scale applied
        %
        %   - axis (string) -
        %     Axis to view the beam from, specified only by 'x', 'y', or, 'z'
        %
        %   - beam_distance (double) -
        %     Distance of the beam from the particle 
        %
        %   - z_sweep (double) -
        %     Range of allowable angles, measured from +z, where the beam
        %     can be translated (rad)
        %
        %   - xy_sweep (double) -
        %     Range of allowable angles, measured from +x towards +y,
        %     where the beam can be translated (rad)
        %
        %   - range (1,2 double) -
        %     Range of values to zoom the image in on
        %
        %   - seed (double) -
        %     Seed to have deterministic randomising of values
            % TODO: Either vary the distance from the particle or beam power
            arguments 
                beam ott.Bsc
                tMatrix ott.Tmatrix
                n_figs single

                % optional parameters
                argv.log_scale logical = false
                argv.beam_distance (1,1) double = 300
                argv.seed (1,1) double = 0
                argv.axis {mustBeText} = 'y'
                argv.range (1,2) double = NaN()
                argv.z_sweep double = pi
                argv.xy_sweep double = pi
            end
            
            % seed randomness
            rng(argv.seed,"twister");

            plots = gobjects(1, n_figs);
            masterFig = figure("Visible", "on");
            master = tiledlayout(masterFig, "flow", ...
                "TileSpacing","compact", ...
                "Padding","compact");
            colormap("gray")

            for i = 1:n_figs

                % randomise beam location
                theta = argv.z_sweep .* rand();
                phi = argv.xy_sweep .* rand();
                transBeam = beam.translateRtp([ argv.beam_distance ; theta ; phi ]);

                % apply T-matrix
                sbeam = tMatrix * transBeam;

                % Figure array
                fig = figure('Visible','off');
                colormap("gray")
                plots(i) = fig;

                % Plot fields
                EColi.plot_scatter(fig, sbeam, argv.axis, argv.range, argv.log_scale);
                EColi.plot_total(fig, sbeam, argv.axis, argv.range, argv.log_scale)

                plotParams = ...
                    sprintf('Parameters: $\\theta = %0.2f\\pi,~ \\phi = %0.2f\\pi$', theta, phi);
                sgtitle(fig, plotParams, 'Interpreter', 'latex')

                % Tiled Figures
                tileFig = nexttile(master);

                sbeam.basis = 'regular';
                if isnan(argv.range)
                    sbeam.visualise('axis', argv.axis, 'log_scale', argv.log_scale);
                else
                    sbeam.visualise('axis', argv.axis, 'range', argv.range, 'log_scale', argv.log_scale);
                end
                caption = sprintf('$\\theta = %0.2f\\pi,~ \\phi = %0.2f\\pi$', ...
                    theta, phi);
                title(tileFig, caption, 'Interpreter', 'latex'); 
                colorbar
            end
        end

        function plots = generate_rand_imgs(beam, wavelength0, n_medium, n_figs, argv)
        % generate_rand_imgs - Generates plots using randomised T-matrices and
        %   beam translations
        %
        % Input Arguments:
        % - beam (Bsc) -
        %   Input beam 
        % 
        % - n_medium (double) -
        %   Refractive index of the medium
        %
        % - wavelength0 (double) -
        %   Initial wavelength0 of the input beam 
        %
        % - n_figs (single) -
        %   Number of desired figures to be generated 
        %
        % - options (name-value pairs) - 
        %   Available options:
        %   
        %   - log_scale (logical) -
        %     Produces output images with a log scale applied
        %
        %   - min_len (double) - 
        %     Minimum length of the E.Coli particle 
        %
        %   - max_len (double) - 
        %     Maximum length of the E.Coli particle
        %
        %   - mean_len (double) - 
        %     Average length of the E.Coli particle
        %
        %   - stdev_len (double) -
        %     Standard deviation for the length of the E.Coli particle
        %
        %   - min_rad (double) - 
        %     Minimum radius of the E.Coli particle 
        %
        %   - max_rad (double) - 
        %     Maximum radius of the E.Coli particle
        %
        %   - mean_rad (double) - 
        %     Average radius of the E.Coli particle
        %
        %   - stdev_rad (double) -
        %     Standard deviation for the radius of the E.Coli particle
        %
        %   - beam_distance (double) -
        %     Distance of the beam from the particle 
        %
        %   - z_sweep (double) -
        %     Range of allowable angles, measured from +z, where the beam
        %     can be translated (rad)
        %
        %   - xy_sweep (double) -
        %     Range of allowable angles, measured from +x towards +y,
        %     where the beam can be translated (rad)
        %
        %   - axis (string) -
        %     Axis to view the beam from, specified only by 'x', 'y', or, 'z'
        %
        %   - range (1,2 double) -
        %     Range of values to zoom the image in on
        %
        %   - seed (double) -
        %     Seed to have deterministic randomising of values

            arguments 
                beam ott.Bsc
                wavelength0 double
                n_medium double
                n_figs single

                % optional parameters
                argv.log_scale logical = false
                argv.min_len double = 1.0
                argv.max_len double = 2.0
                argv.mean_len double = 4.0
                argv.stdev_len double = 0.08
                argv.min_rad double = 0.35
                argv.max_rad double = 1.0
                argv.mean_rad double = 0.89
                argv.stdev_rad double = 0.08
                argv.beam_distance (1,1) double = 300
                argv.axis {mustBeText} = 'y'
                argv.range (1,2) double = NaN()
                argv.seed double = 0
                argv.z_sweep double = pi
                argv.xy_sweep double = pi
            end

            % seed randomiser
            rng(argv.seed, "twister")

            % normal distirbution of len and rad
            len = argv.mean_len + argv.stdev_len .* randn(1, n_figs);
            rad = argv.mean_rad + argv.stdev_rad .* randn(1, n_figs);

            % allocate plots
            plots = gobjects(1, n_figs);
            masterFig = figure("Visible", "on");
            master = tiledlayout(masterFig, "flow", ...
                "TileSpacing","compact", ...
                "Padding","compact");
            colormap("gray")

            for i = 1 : n_figs

                % generate T-matrix
                T = EColi.T_matrix(wavelength0, n_medium, ...
                    "length", len(i), ...
                    "radius", rad(i) ...
                );

                % randomise beam location
                theta = argv.z_sweep .* rand();
                phi = argv.xy_sweep .* rand();
                transBeam = beam.translateRtp([ argv.beam_distance ; theta ; phi ]);

                % apply T-matrix
                sbeam = T * transBeam;

                % Figure array
                fig = figure('Visible','off');
                colormap("gray")
                plots(i) = fig;

                % Plot fields
                EColi.plot_scatter(fig, sbeam, argv.axis, argv.range, argv.log_scale);
                EColi.plot_total(fig, sbeam, argv.axis, argv.range, argv.log_scale)

                plotParams = ...
                    sprintf('Parameters: $\\ell = %0.2f,~ r = %0.2f,~ \\theta = %0.2f\\pi,~ \\phi = %0.2f\\pi$', ...
                        len(i), rad(i), theta, phi);
                sgtitle(fig, plotParams, 'Interpreter', 'latex')

                % Tiled Figures
                tileFig = nexttile(master);

                sbeam.basis = 'regular';
                if isnan(argv.range)
                    sbeam.visualise('axis', argv.axis, 'log_scale', argv.log_scale);
                else
                    sbeam.visualise('axis', argv.axis, 'range', argv.range, 'log_scale', argv.log_scale);
                end

                caption = sprintf(...
                    '$\\theta = %0.2f\\pi,~ \\phi = %0.2f\\pi,~ \\ell = %0.2f,~ r = %0.2f$', ...
                    theta, phi, len(i), rad(i) ...
                );
                title(tileFig, caption, 'Interpreter', 'latex'); 
                colorbar
            end
        end

        function test_Tmatrix(tMatrix, tol)
        %TEST_TMATRIX - runs a series of tests on a given T-matrix
        % Input Arguments:
        % - tMatrix (T_matrix) -
        %   Input T-matrix to be tested on
        %
        % - tol (double) -
        %   Tolerance for power testing
            EColi.testPower(tMatrix, tol)
            EColi.testConvert(tMatrix)
            EColi.testResizing(tMatrix)
            EColi.testRealImageFunctions(tMatrix)
        end
    end

    methods (Access = private, Static)
        function plot_scatter(fig, sbeam, axis, range, log)
            sbeam.basis = 'regular';
            subplot(1,2,1, 'Parent', fig);

            if isnan(range)
                sbeam.visualise('axis', axis, 'log_scale', log);
            else
                sbeam.visualise('axis', axis, 'range', range, 'log_scale', log);
            end

            title('Scattered Field')
            colorbar
        end

        function plot_total(fig, sbeam, axis, range, log)
            sbeam.basis = 'outgoing';

            tbeam = sbeam.totalField(sbeam);
            tbeam.basis = 'regular';

            subplot(1,2,2, 'Parent', fig);
            if isnan(range)
                tbeam.visualise('axis', axis, 'log_scale', log);
            else
                tbeam.visualise('axis', axis, 'range', range, 'log_scale', log);
            end

            title('Total Field')
            colorbar
        end

        function [len, rad] = randomise_size(min_len, max_len, min_rad, max_rad)
            len = min_len + (max_len - min_len) .* rand();
            rad = min_rad + (max_rad - min_rad) .* rand();
        end

        function [rho, z] = parameterise_shape(n_pts, R, H)
            arguments 
                n_pts single 
                R double
                H = double
            end

            capPts = n_pts / 2;

            angle = linspace(0, pi/2, capPts);
            rhoLeft = R * sin(angle);
            zLeft = (H / 2) + R * cos(angle);

            angle = linspace(pi/2, pi, capPts);
            rhoRight = R * sin(angle);
            zRight = -(H / 2) + R * cos(angle);

            rho = [rhoLeft, rhoRight];
            z = [zLeft, zRight];
        end

        function testPower(T, tol)
            try
            EColi.checkWithTol(abs(1.0 - sum(abs(T.total.data).^2, 2)), ... 
                tol, ...
                'Power not conserved');
            catch
                colSum = sum(abs(T.total.data).^2, 1);   % column-wise sum
                deviation = abs(1.0 - colSum);
                offending = deviation > tol;

                fprintf('Largest column sum deviation:\n')
                if any(offending)
                    disp(max(deviation(offending)))
                else
                    disp(0)
                end

                fprintf('Smallest column sum deviation:\n')
                if any(offending)
                    disp(min(deviation(offending)))
                else
                    disp(0)
                end

                fprintf('Number of offending columns: %d\n', sum(offending))
            end

            tTotal = T.total;
            tTotal.Nmax = tTotal.Nmax + 5;
            EColi.checkWithTol(abs(1.0 - sum(abs(tTotal.data).^2, 2)), ...
                tol, ...
                'Resizing T-total does not conserve power');

            tScat = T.scattered;
            tScat.Nmax = tScat.Nmax + 5;
            EColi.checkWithTol(abs(1.0 - sum(abs(tScat.total.data).^2, 2)), ...
                tol, ...
                'Resizing T-total does not conserve power');
        end

        function checkWithTol(val, tol, msg)
            assert(all(val < tol), msg);
            % if (val < tol)
            %     bool = true;
            % else
            %     bool = false; 
            %     fprintf(msg)
            % end
        end

        function testConvert(T)
            total = T.total;
            scat = T.scattered;

            assert(strcmpi(total.type, 'total'), 'Conversion to total');
            assert(strcmpi(scat.type, 'scattered'), ...
                'Conversion to scattered');

            % Check conversions back
            assert(strcmpi(scat.total.type, 'total'), ...
                'Conversion back to total');
            assert(strcmpi(total.scattered.type, 'scattered'), ...
                'Conversion back to scattered');
        end

        function testResizing(T)
            Tnew1 = T;

            Tnew1.Nmax = Tnew1.Nmax + 5;
            assert(all(Tnew1.Nmax == T.Nmax + 5), ...
                'Failed to increase Nmax with vector size');
            assert(all(size(Tnew1.data) > size(T.data)), ...
                'Tmatrix size not actually increased (vector input)');

            Tnew2 = Tnew1;
            Tnew2.Nmax = T.Nmax;
            assert(all(Tnew2.Nmax == T.Nmax), ...
                'Failed to decrease Nmax with vector size');
            assert(all(size(Tnew2.data) == size(T.data)), ...
                'Tmatrix size not actually decreased (vector input)');

            Tnew1 = T;
            Tnew1.Nmax = T.Nmax(1) + 5;
            assert(all(Tnew1.Nmax == T.Nmax + 5), ...
                'Faild to increase Nmax (scalar input)');
            assert(all(size(Tnew1.data) > size(T.data)), ...
                'Tmatrix size not actually increased (scalar input)');

            Tnew1 = T;
            Tnew1.Nmax = [T.Nmax(1) + 5, T.Nmax(2)];
            assert(all(Tnew1.Nmax == [T.Nmax(1) + 5, T.Nmax(2)]), ...
                'Failed to increase Nmax (uneven input)');
            assert(size(Tnew1.data, 1) > size(T.data, 1) ...
                && size(Tnew1.data, 2) == size(T.data, 2), ...
                'Tmatrix size not increased correctly (uneven input)');

            Tnew1 = T;
            Tnew1.Nmax(1) = T.Nmax(1) + 5;
            assert(all(Tnew1.Nmax == [T.Nmax(1) + 5, T.Nmax(2)]), ...
                'Faild to increase Nmax (index input)');
            assert(size(Tnew1.data, 1) > size(T.data, 1) ...
                && size(Tnew1.data, 2) == size(T.data, 2), ...
                'Tmatrix size not increased correctly (index input)');
        end

        function testRealImageFunctions(T)
            T2 = real(T);
            T3 = imag(T);
            T4 = T.real();
            T5 = T.imag();
        end
    end
end
