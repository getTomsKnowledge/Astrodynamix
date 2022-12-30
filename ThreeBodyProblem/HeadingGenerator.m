classdef HeadingGenerator
    properties
        courseId
        assignmentNum
        description
    end
    methods
        function newObj = HeadingGenerator(id, num, descrip)

            newObj.courseId = id;
            newObj.assignmentNum = num;
            newObj.description = descrip;

        end       
        %void printHeading()
        function print_heading(obj)

            formatting = PrintFormatting();
            fprintf('%s\n%s\n\n', ...
                formatting.lineBreak, formatting.lineBreak);
            fprintf('   Name:        %s\n', formatting.myName);
            fprintf('   Course:      %s\n', obj.courseId);
            fprintf('   Assignment:  %s\n', obj.assignmentNum);
            fprintf('   Description: %s\n\n', obj.description);
            fprintf('%s\n%s\n\n', ...
                formatting.lineBreak, formatting.lineBreak);

        end
    end
    
end