%  
classdef parseInputArgsObject
    methods(Access = public, Static, Hidden = true)
        % parse
        function p = parseInputArgs(varargin)
            p = inputParser;
            addParamValue(p,'region',[]);
            addParamValue(p,'integrand',[]);
            addParamValue(p,'integrandInput',[]);
            addParamValue(p,'outerRegion',[]);
            addParamValue(p,'innerRegion',[]);
            addParamValue(p,'region1',[]);
            addParamValue(p,'region2',[]);
            addParamValue(p,'name',[]);
            addParamValue(p,'p',[]);
            addParamValue(p,'e',[]);
            addParamValue(p,'t',[]);
            addParamValue(p,'u',[]);
            addParamValue(p,'time',[]);
            addParamValue(p,'geometry',[]);
            addParamValue(p,'fConstantValue',[]);
            addParamValue(p,'fiFunction',[]);
            addParamValue(p,'cConstantValue',[]);
            addParamValue(p,'aConstantValue',[]);
            addParamValue(p,'dConstantValue',[]);
            addParamValue(p,'cijFunction',[]);
            addParamValue(p,'aijFunction',[]);
            addParamValue(p,'dijFunction',[]);
            addParamValue(p,'triangleIndividualCoordinatesForm',[]);
            addParamValue(p,'symbolicEquationFunction',[]);
            addParamValue(p,'symbolicFunction',[]);
            addParamValue(p,'forceGen',true);
            addParamValue(p,'parameters',[]);
            addParamValue(p,'xyutFunction',[]);
            addParamValue(p,'dirichlet',[]);
            addParamValue(p,'neumann',[]);
            addParamValue(p,'type',[]);
            addParamValue(p,'nodes',[]);
            addParamValue(p,'value',[]);
            addParamValue(p,'dimension',[]);
            addParamValue(p,'center',[]);
            addParamValue(p,'radius',[]);
            addParamValue(p,'image',[]);
            addParamValue(p,'minimumPoints',[]);
            addParamValue(p,'clockwise',[]);
            addParamValue(p,'leaveOpen',[]);
            addParamValue(p,'edgeStartPoint',[]);
            addParamValue(p,'edgeEndPoint',[]);           
            addParamValue(p,'points',[]);
            addParamValue(p,'leftRegion',[]);
            addParamValue(p,'rightRegion',[]);
            addParamValue(p,'leftRegionIsExterior',[]);
            addParamValue(p,'interiorPt',[]);
            addParamValue(p,'periodicX',[]);
            addParamValue(p,'periodicY',[]);
            
            addParamValue(p,'showMesh',false);
            addParamValue(p,'showNodes',false);
            addParamValue(p,'displayCoefficients',false);
            addParamValue(p,'showBoundary',false);
            addParamValue(p,'showBoundaryClear',false);
            addParamValue(p,'doNotSum',false);
            addParamValue(p,'numRefineMeshSteps',0);

            p.KeepUnmatched = true;
            parse(p,varargin{:});
        end
    end
end