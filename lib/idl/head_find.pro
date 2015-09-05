function head_find, fname, keyword, ext=ext

    if not n_params() eq 2 then begin
        print, 'value = head_find(fname, keyword, ext=ext)'
        print, 'Find value of "keyword" in header of P-Smac2'
        print, 'generated fitsfile'
    end

    if not keyword_set(ext) then $
        ext = 0

    ; Some keywords are not in the comments
    if strcmp(keyword, 'NAXIS1')                or $
       strcmp(keyword, 'NAXIS2')                or $
       strcmp(keyword, 'PROG')                  or $
       strcmp(keyword, 'CITE')                  or $
       strcmp(keyword, 'VERSION')               or $
       strcmp(keyword, 'NTASK')                 or $
       strcmp(keyword, 'COMPILETIMEOPTIONS')    or $ 
       strcmp(keyword, 'MAKEFILE')              or $
       strcmp(keyword, 'Unit')                  or $
       strcmp(keyword, 'Module Name')           or $
       strcmp(keyword, 'NAME')                  or $
       strcmp(keyword, 'Description')           or $
       strcmp(keyword, 'DESCR')                 then begin
        substridx = 1
    end else begin ;parameters are in the comments
        substridx = 0
    end

    img = readfits(fname, header,/silent, ext=ext)
    img = 0

    ; find first occurance
    for i=0,n_elements(header)-1 do begin
        
        str_found = strpos(header[i], keyword)
        
        if str_found ne -1 then $
            break
    end

    ; catch not found
    if str_found eq -1 then begin
        print, 'Keyword  not found : ', keyword
        return, -1
    end

    ; extract value or comment
    substrn = strsplit(header[i],'=',/extract)  
    type = substrn[0]
    ssubstrn = strsplit(substrn[1], '/',/extract)
    result = ssubstrn[substridx]

    ; treat special cases
    if strcmp(keyword, 'Unit')  then $
        for i=2,n_elements(ssubstrn)-1 do $
           result += '/'+ssubstrn[i] 

    ; treat type
    if strcmp(type, 'INT')      then $
        result = long(result)

    if strcmp(type, 'DOUBLE')   then $
        result = double(result)    

    return, result
end
