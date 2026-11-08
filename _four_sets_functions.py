def first_stage_mask(sx, cy, outward, start_offset, start_width, end_width, length=40.0):
    """Smooth full-width release envelope, with exact widths at both ends."""
    import klayout.db as kdb
    import numpy as np
    require(0 <= start_offset < length, 'Invalid first-stage start offset')
    require(start_width >= end_width > 0, 'Mask must narrow monotonically')
    s = np.linspace(start_offset, length, 401)
    t = (s - start_offset) / (length - start_offset)
    f = t**3 * (10 - 15*t + 6*t*t)
    w = start_width + (end_width-start_width)*f
    upper = [(sx+outward*float(a), cy+float(b)/2) for a,b in zip(s,w)]
    lower = [(sx+outward*float(a), cy-float(b)/2) for a,b in zip(s[::-1],w[::-1])]
    return kdb.Region(kdb.Polygon([kdb.Point(round(x/DBU),round(y/DBU)) for x,y in upper+lower]))


def finalize_coupler_windows(layout, entries):
    """Apply the independent entrance widths, start positions and margins.

    L2 starts 0.5 um outward from L10, preserving the original delayed-copy
    interface. L14 starts with L10. Entrance widths are TOTAL widths, not
    extra clearances: L2=15 um and L14=8 um. Both downstream margins=2.5 um.
    Only the first 40 um window envelope is replaced; downstream mask
    candidates still come from the original mask generator.
    """
    import klayout.db as kdb
    policies = {}
    require(COUPLER_END_WINDOWS_ONLY, 'This comparison requires terminal-only masks')
    for e in entries:
        cell = layout.cell(e['case_id'])
        rr = regions(layout,cell)
        ml,mc = read_gds(WORK/e['release_window_file'])
        candidates = regions(ml,mc)
        result = {2:kdb.Region(),14:kdb.Region()}
        ends=[]
        for side,port in (('input',e['input_port_local_um']),('output',e['output_port_local_um'])):
            px,cy = map(float,port)
            sign = -1 if px<0 else 1
            search = box((-20000,cy/8,px,20000) if sign<0 else (px,-20000,20000,cy/8))
            taper = (rr[(10,0)] & search).merged()
            if taper.is_empty():
                require(side=='output' and e['mode']=='straight_output','Unexpected absent taper')
                continue
            sx = (taper.bbox().right if sign<0 else taper.bbox().left)*DBU
            end = dict(side=side,center_y_um=cy,l10_start_local_x_um=sx)
            for layer,offset,start_width,margin in (
                (2,LAYER2_START_OFFSET,LAYER2_START_WIDTH,LAYER2_WINDOW_MARGIN),
                (14,LAYER14_START_OFFSET,LAYER14_START_WIDTH,LAYER14_WINDOW_MARGIN)):
                start = sx+sign*offset
                keep = box((-20000,cy/8,start,20000) if sign<0 else (start,-20000,20000,cy/8))
                retained = (rr[(layer,0)] & candidates[(layer,0)] & keep).merged()
                endpoint = sx+sign*40.0
                replace = box((min(sx,endpoint),cy-25,max(sx,endpoint),cy+25))
                shaped = first_stage_mask(sx,cy,sign,offset,start_width,1.0+2*margin)
                retained = ((retained-replace) | shaped).merged()
                require((retained-keep).is_empty(),'Mask extends before its entrance')
                result[layer] |= retained
                end[f'layer{layer}_start_local_x_um']=start
                end[f'layer{layer}_start_global_x_um']=start+e['origin_um'][0]
                end[f'layer{layer}_start_width_um']=start_width
                end[f'layer{layer}_offset_um']=offset
            ends.append(end)
        require(len(ends)==(1 if e['mode']=='straight_output' else 2),'Wrong number of couplers')
        for n,r in result.items():
            cell.shapes(layout.layer(n,0)).clear()
            cell.shapes(layout.layer(n,0)).insert(r.merged())
        policies[e['case_id']]=dict(ends=ends,layer2_margin_um=LAYER2_WINDOW_MARGIN,
            layer14_margin_um=LAYER14_WINDOW_MARGIN)
        e.update(final_pnc_present=False,final_release_window_count=len(ends),
                 layer2_margin_um=LAYER2_WINDOW_MARGIN,layer14_margin_um=LAYER14_WINDOW_MARGIN,
                 layer2_start_offset_um=LAYER2_START_OFFSET,layer14_start_offset_um=LAYER14_START_OFFSET,
                 layer2_start_width_um=LAYER2_START_WIDTH,layer14_start_width_um=LAYER14_START_WIDTH)
    return policies


def make_single_straight(source_regions, width, single_coupler):
    """Copy complete approved terminal assemblies and replace only routing.

    New ports are (-3000,0) and (+3000,0). Move the complete upper-left end
    down 200 um and lower-right end up 200 um; no optical taper is redrawn.
    The one-coupler version retains the 2000 um taper-free output extension.
    """
    import klayout.db as kdb
    left = box((-20000,25,-STRAIGHT_LENGTH/2,20000))
    right = box((STRAIGHT_LENGTH/2,-20000,20000,-25))
    down = kdb.Trans(0,round(-200/DBU))
    up = kdb.Trans(0,round(200/DBU))
    out={}
    for key,r in source_regions.items():
        if key[0] in REMOVED_OUTPUT_LAYERS:
            continue
        out[key]=((r&left).transformed(down) | (r&right).transformed(up)).merged()
    ridge=box((-STRAIGHT_LENGTH/2,-width/2,STRAIGHT_LENGTH/2,width/2))
    out[(15,0)]=(out.get((15,0),kdb.Region())|ridge).merged()
    return {key:r for key,r in out.items() if not r.is_empty()}


def combine():
    """Place four sets of three widths; preserve the first six device origins."""
    import copy
    import gdsfactory as gf
    import klayout.db as kdb
    gf.gpdk.PDK.activate()
    base=load_base()
    ly=kdb.Layout();ly.dbu=DBU
    top=ly.create_cell('D1_L6000_STRAIGHT_OUT_COMPARISON')
    entries=[];checks=[];bounds=kdb.Box();texts=[]

    def place(ident,rr,entry):
        nonlocal bounds
        row=len(entries);y=FIRST_ROW_Y+row*ROW_PITCH
        x=entry['origin_um'][0]
        cell=insert_regions(ly,ident,rr)
        tr=kdb.Trans(round(x/DBU),round(y/DBU))
        top.insert(kdb.CellInstArray(cell.cell_index(),tr))
        entry.update(case_id=ident,origin_um=[tr.disp.x*DBU,tr.disp.y*DBU])
        for (n,d),r in rr.items():
            if n in FAB_LAYERS:bounds+=r.transformed(tr).bbox()
        text=gf.components.text(text=entry['label'],size=LABEL_SIZE,layer=(2,0))
        textreg=kdb.Region(text.kdb_cell.begin_shapes_rec(gf.get_layer((2,0)))).merged()
        labelcell=insert_regions(ly,'LABEL_'+ident,{(2,0):textreg})
        lt=kdb.Trans(round((LEFT_BEND_X+200)/DBU),round((y+250)/DBU))
        top.insert(kdb.CellInstArray(labelcell.cell_index(),lt))
        texts.append(textreg.transformed(lt));bounds+=texts[-1].bbox()
        entry['label_origin_um']=[lt.disp.x*DBU,lt.disp.y*DBU]
        entries.append(entry)

    # A/B: existing meander designs, retaining every original end and origin.
    modes=([False] if INCLUDE_REFERENCE_DEVICES else [])+[True]
    for one in modes:
        for width in WG_WIDTHS:
            meta=json.loads((WORK/(base.case_id(width)+'.json')).read_text())
            ident=meta['variant_id'] if one else meta['reference_id']
            src,sc=read_gds(WORK/(ident+'.gds'))
            rr=regions(src,sc)
            label=f'D1_W{width:.2f}_T80_P10'+(f'_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '')
            place(ident,rr,dict(mode='straight_output' if one else 'original_two_tapers',
                topology='meander',set_name='meander_one_coupler' if one else 'meander_two_couplers',
                width_um=width,origin_um=[LEFT_BEND_X-meta['left_pnc_x_um'],0],label=label,
                straight_section_count=3,straight_length_um=STRAIGHT_LENGTH,
                bend_count=2,added_length_um=OUTPUT_EXTENSION_LENGTH if one else 0.,
                path_length_excluding_tapers_um=meta['meander_length_um']+(OUTPUT_EXTENSION_LENGTH if one else 0.),
                release_window_file=meta['release_window_file'],
                input_port_local_um=meta['input_port_um'],output_port_local_um=meta['output_port_um']))
    policies=finalize_coupler_windows(ly,entries)

    # C/D: one straight only; derive complete couplers from the updated A/B cells.
    for one in (False,True):
        for width in WG_WIDTHS:
            source_id=base.case_id(width)+(f'_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '')
            require(ly.cell(source_id) is not None,'Enable reference devices for the straight two-coupler set')
            source_entry=next(e for e in entries if e['case_id']==source_id)
            rr=make_single_straight(regions(ly,ly.cell(source_id)),width,one)
            ident=f'D1_SINGLE_L{STRAIGHT_LENGTH:g}_W{round(width*1000)}'+(f'_1T_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '_2T')
            entry=copy.deepcopy(source_entry)
            entry.update(mode='single_straight_one_coupler' if one else 'single_straight_two_couplers',
                topology='single_straight',set_name='single_straight_one_coupler' if one else 'single_straight_two_couplers',
                straight_section_count=1,bend_count=0,
                path_length_excluding_tapers_um=STRAIGHT_LENGTH+(OUTPUT_EXTENSION_LENGTH if one else 0.),
                input_port_local_um=[-STRAIGHT_LENGTH/2,0.],output_port_local_um=[STRAIGHT_LENGTH/2,0.],
                label=f'D1_S{STRAIGHT_LENGTH:g}_W{width:.2f}_T80'+(f'_1T_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '_2T'),
                copied_from=source_id)
            place(ident,rr,entry)
            policies[ident]=copy.deepcopy(policies[source_id])
            for end in policies[ident]['ends']:end['center_y_um']=0.

    # Retain original outline width; extend height to contain the added rows.
    xmin=math.floor((bounds.left*DBU-FRAME_MARGIN)/100)*100
    xmax=math.ceil((bounds.right*DBU+FRAME_MARGIN)/100)*100
    ymin=math.floor((bounds.bottom*DBU-FRAME_MARGIN)/100)*100
    ymax=math.ceil((bounds.top*DBU+FRAME_MARGIN)/100)*100
    model=base.configured_design(WG_WIDTHS[0],dict(CHIP_WIDTH=xmax-xmin,CHIP_MIN_HEIGHT=max(MIN_FRAME_HEIGHT,ymax-ymin)))
    actual_height=max(MIN_FRAME_HEIGHT,ymax-ymin,(len(entries)-1)*model.DEVICE_Y_PITCH+model.CHIP_Y_ALLOWANCE)
    ymax=ymin+actual_height
    frame=gf.Component();model.add_chip_outline_and_markers(frame,len(entries));frame.flatten()
    fc=insert_regions(ly,'FRAME_AND_MARKERS',regions(frame.kcl.layout,frame.kdb_cell))
    ft=kdb.Trans(round((xmin+xmax)/2/DBU),round((ymin+ymax)/2/DBU))
    top.insert(kdb.CellInstArray(fc.cell_index(),ft))
    layer_policy=finalize_output_layers(ly,fc)
    snapshots={e['case_id']:regions(ly,ly.cell(e['case_id'])) for e in entries}
    placed=[{n:r.transformed(kdb.Trans(round(e['origin_um'][0]/DBU),round(e['origin_um'][1]/DBU)))
             for (n,d),r in snapshots[e['case_id']].items()} for e in entries]
    for i,a in enumerate(placed):
        for j,b in enumerate(placed[i+1:],i+1):
            area=sum((a.get(n,kdb.Region())&b.get(n,kdb.Region())).area()*DBU**2 for n in (2,10,13,14,15))
            check(checks,f'rows_{i}_{j}_no_overlap',area==0,area_um2=area)
    for i,t in enumerate(texts):
        check(checks,f'label{i}_clear',all((t&r).is_empty() for p in placed for r in p.values()))
        check(checks,f'label{i}_unique',all((t&s).is_empty() for s in texts[i+1:]))
    fr=regions(ly,fc)
    for n in (10,15):
        marks=fr[(n,0)].transformed(ft)
        check(checks,f'L{n}_markers_clear',all((marks&r).is_empty() for p in placed for r in p.values()))
    write_layout(ly,OUTPUT)
    final,finaltop=read_gds(OUTPUT)
    insts={i.cell.name:i for i in finaltop.each_inst()}
    check(checks,'twelve_devices_and_labels',len(entries)==12 and len(insts)==25)
    check(checks,'excluded_layers_absent',not ({final.get_info(i).layer for i in final.layer_indexes()}&set(REMOVED_OUTPUT_LAYERS)))
    for e in entries:
        rr=regions(final,insts[e['case_id']].cell)
        check(checks,e['case_id']+'_saved_layers',set(rr)==set(snapshots[e['case_id']]))
        for key,r in snapshots[e['case_id']].items():
            check(checks,e['case_id']+f'_L{key[0]}_export_exact',(r^rr[key]).is_empty())
    workers=[json.loads((WORK/(base.case_id(w)+'.json')).read_text()) for w in WG_WIDTHS]
    report=dict(status='PASS',devices=entries,device_count=len(entries),set_count=4,
        width_sweep_um=list(WG_WIDTHS),meander_straight_length_um=STRAIGHT_LENGTH,
        single_straight_length_um=STRAIGHT_LENGTH,extension_length_um=OUTPUT_EXTENSION_LENGTH,
        row_pitch_um=ROW_PITCH,left_bend_x_um=LEFT_BEND_X,frame_bounds_um=[xmin,ymin,xmax,ymax],
        layer2_start_width_um=LAYER2_START_WIDTH,layer14_start_width_um=LAYER14_START_WIDTH,
        layer2_window_margin_um=LAYER2_WINDOW_MARGIN,layer14_window_margin_um=LAYER14_WINDOW_MARGIN,
        layer2_start_offset_um=LAYER2_START_OFFSET,layer14_start_offset_um=LAYER14_START_OFFSET,
        final_layer_policy=layer_policy,coupler_layer_policy=policies,checks=checks,workers=workers,
        check_count=len(checks)+sum(len(w['checks']) for w in workers),
        original_gds_sha256=hashlib.sha256(BASE_GDS.read_bytes()).hexdigest(),
        gds_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        optical_or_mechanical_simulation_performed=False)
    OUTPUT.with_suffix('.json').write_text(json.dumps(report,indent=2))
    with OUTPUT.with_suffix('.csv').open('w',newline='') as fp:
        fields=['case_id','set_name','width_um','origin_um','straight_section_count','bend_count','added_length_um','path_length_excluding_tapers_um','label']
        writer=csv.DictWriter(fp,fieldnames=fields);writer.writeheader()
        writer.writerows({k:e[k] for k in fields} for e in entries)
    print('PASS',report['check_count'],'construction/export checks;',OUTPUT,flush=True)
