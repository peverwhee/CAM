#! /bin/bash

# Tool to produce baselines useful for validating the move to NUOPC tests
# Convert tests to version including driver (mct, nuopc) name
# A new directory is created for most tests from a baseline and MCT baselines
#   are modified to look as if it was a NUOPC baseline.
# In addition, some tests also have 'copies' which is just a softlink to the
#   original baseline

MCT_CP=""
MCT_CP="${MCT_CP} ERC_D_Ln9_Vmct.f10_f10_mg37.FHS94.izumi_nag.cam-idphys"
MCT_CP="${MCT_CP} ERP_Ln9_Vmct.ne5pg3_ne5pg3_mg37.QPC6.izumi_nag.cam-outfrq9s_clubbmf"
MCT_CP="${MCT_CP} "
MCT_CP="${MCT_CP} "
MCT_CP="${MCT_CP} "
MCT_CP="${MCT_CP} "

MESH_MAP=""
MESH_MAP="${MESH_MAP} ne5pg3:ne5pg3_ESMFmesh_cdf5_c20210118.nc"
MESH_MAP="${MESH_MAP} ne16:ne16np4_scrip_171002_ESMFmesh.nc"
MESH_MAP="${MESH_MAP} ne16pg3:ne16pg3_ESMFmesh_cdf5_c20211018.nc"
MESH_MAP="${MESH_MAP} ne5:ne5np4_esmf_mesh_c210121.nc"
MESH_MAP="${MESH_MAP} ne5pg3:ne5pg3_ESMFmesh_cdf5_c20210118.nc"
MESH_MAP="${MESH_MAP} T42:T42_ESMFmesh_c20200629.nc"
MESH_MAP="${MESH_MAP} T5:T5_ESMFmesh_cdf5_c20210923.nc"
MESH_MAP="${MESH_MAP} ne5pg4:ne5pg4_ESMFmesh_cdf5_c20210118.nc"
MESH_MAP="${MESH_MAP} ne5pg2:ne5pg2_ESMFmesh_cdf5_c20211018.nc"

if [ "${HOST:0:5}" == "izumi" ]; then
  machname="Izumi"
  BL_ROOT="/fs/cgd/csm/models/atm/cam/pretag_bl"
  nuopc_store="/home/courtneyp/baselines/output/nuopc_store"
  mesh_base="/fs/cgd/csm/inputdata/share/meshes/"
  if [ -z "${dest_root}" ]; then
    dest_root="${BL_ROOT}"
  fi
elif [ "${HOST:0:8}" == "cheyenne" ]; then
  machname="Cheyenne"
  BL_ROOT="/glade/p/cesm/amwg/cesm_baselines"
  mesh_base="/glade/p/cesm/cseg/inputdata/share/meshes/"
  nuopc_store="??"
  if [ -z "${dest_root}" ]; then
    dest_root="${BL_ROOT}"
  fi
else
  perr "ERROR: Unsupported host, \"${HOST}\"." 3
fi

perr() {
  if [ $# -ge 1 ]; then
    echo "${1}"
  fi
  if [ $# -ge 2 ]; then
    exit $2
  else
    exit 1
  fi
}

softlink_file() {
  # Args: <filename>, <dir_orig>, <dir_dest>
  local tfile="${1}"
  local src_dir="${2}"
  local dst_dir="${3}"
  ln -s ${src_dir}/${tfile} ${dst_dir}/${tfile}
}

copy_dir() {
  # Args: <subdir_name>, <dir_orig>, <dir_dest>
  local dir_name="${1}"
  local dir_orig="${2}/${dir_name}"
  local dir_dest="${3}/${dir_name}"
  if [ ! -d "${dir_dest}" ]; then
    mkdir ${dir_dest}
  fi
  for file in `cd ${dir_orig}; ls`; do
    cp ${dir_orig}/${file} ${dir_dest}/${file}
  done
}

write_med_modelio() {
  # Write med_modelio.nml to $1
  cat > ${1} << EOF
&pio_inparm
  pio_netcdf_format = "64bit_offset"
  pio_numiotasks = -99
  pio_rearranger = 2
  pio_root = 1
  pio_stride = 48
  pio_typename = "netcdf"
/

EOF
}

write_drv_in() {
  # Write drv_in to $1
  cat > ${1} << EOF
&debug_inparm
  create_esmf_pet_files = .false.
/
&papi_inparm
  papi_ctr1_str = "PAPI_FP_OPS"
  papi_ctr2_str = "PAPI_NO_CTR"
  papi_ctr3_str = "PAPI_NO_CTR"
  papi_ctr4_str = "PAPI_NO_CTR"
/
&pio_default_inparm
  pio_async_interface = .false.
  pio_blocksize = -1
  pio_buffer_size_limit = -1
  pio_debug_level = 0
  pio_rearr_comm_enable_hs_comp2io = .true.
  pio_rearr_comm_enable_hs_io2comp = .false.
  pio_rearr_comm_enable_isend_comp2io = .false.
  pio_rearr_comm_enable_isend_io2comp = .true.
  pio_rearr_comm_fcd = "2denable"
  pio_rearr_comm_max_pend_req_comp2io = -2
  pio_rearr_comm_max_pend_req_io2comp = 64
  pio_rearr_comm_type = "p2p"
/
&prof_inparm
  profile_add_detail = .false.
  profile_barrier = .false.
  profile_depth_limit = 12
  profile_detail_limit = 2
  profile_disable = .false.
  profile_global_stats = .true.
  profile_outpe_num = 1
  profile_outpe_stride = 0
  profile_ovhd_measurement = .false.
  profile_papi_enable = .false.
  profile_single_file = .false.
  profile_timer = 4
/
EOF
}

add_meshfile_to_docn_in() {
  local docn_file="${2}/CaseDocs/docn_in"
  local temp_file="${2}/CaseDocs/temp.txt"
  local test_name="$(basename ${2})"
  local count=1
  local mesh=""
  local meshfile=""
  local found=0
  local mesh_trimmed=""
  local mesh_available=0
  # Get mesh from test_name
  old_ifs="$IFS"
  IFS="."
  for i in ${test_name}; do
    if [ ${count} -eq 2 ]; then
      IFS="_"
      for j in ${i}; do
        mesh="${i}"
      done
      break
    fi
    ((count=count+1))
  done
  IFS="_"
  for i in ${mesh}; do
    mesh_trimmed="${i}"
    break
  done
  IFS=" "
  for mesh_map_i in ${MESH_MAP}; do
    IFS=":"
    for i in ${mesh_map_i}; do
      if [ ${found} -eq 1 ]; then
        meshfile="`echo ${i} | xargs`"
        found=0
        break
      fi
      trimmed_i="`echo ${i} | xargs`"
      if [ "${trimmed_i}" == "${mesh_trimmed}" ]; then
        found=1
        mesh_available=1
      else
        found=0
      fi
    done
    IFS=" "
  done
  if [ ${mesh_available} -eq 1 ] && [ -f ${docn_file} ]; then
    # Write docn_n
#    cat "${docn_file}" | head -2 > "${temp_file}"
#    echo "  model_maskfile=\""${mesh_base}${meshfile}"\"" >> ${temp_file}
#    echo "  model_meshfile=\""${mesh_base}${meshfile}"\"" >> ${temp_file}
#    cat "${docn_file}" | tail -5 >> "${temp_file}"
    cat > ${docn_file} << EOF
&docn_nml
  datamode = "sst_aquap3"
  model_maskfile = "${mesh_base}${meshfile}"
  model_meshfile = "${mesh_base}${meshfile}"
  nx_global = 1352
  ny_global = 1
  restfilm = "null"
  sst_constant_value = -1.0
/

EOF
  fi
  IFS="$old_ifs"
}

remove_unnecessary_nl_entries() {
  local test_file="${2}/CaseDocs/${1}"
  local temp_file="${2}/CaseDocs/temp.nml"
  cat "${test_file}" | tail -8 > "${temp_file}"
  mv ${temp_file} "${test_file}"
}

clone_bl_test_dir() {
  # Args: <orig_bl> <driver_bl> <bl_test_dir> <driver>
  local test_dir="${3}"
  local src_bl_tdir="${1}/${test_dir}"
  local driver="${4}"
#  local fmod="s/[.]/_V${driver}./"
  local fmod="s/Vmct/V${driver}/"
  # We need the new driver name of the test directory
  nuopc_dir="`echo ${test_dir} | sed -e ${fmod}`"
  local dest_bl_tdir="${2}/${nuopc_dir}"
  local tfile
  # Create the destination test directory
  if [ ! -d "${dest_bl_tdir}" ]; then
      mkdir "${dest_bl_tdir}"
    fi
  # Create soft links for most files
  # Looking for *.nc user* TestStatus *.log* ?
  for tfile in `cd ${src_bl_tdir}; ls`; do
    if [ -f "${src_bl_tdir}/${tfile}" ] && [[ ${tfile} != "cpl.hi"* ]]; then
      softlink_file "${tfile}" "${src_bl_tdir}" "${dest_bl_tdir}"
    fi
  done
  # Copy contents of CaseDocs
  copy_dir "CaseDocs" "${src_bl_tdir}" "${dest_bl_tdir}"
  for tfile in `cd ${dest_bl_tdir}/CaseDocs; ls *.nml`; do
    #if [ -f "${src_bl_tdir}/CaseDocs/${tfile}" ]; then
#      remove_unnecessary_nl_entries "${tfile}" "${dest_bl_tdir}"
      write_med_modelio "${dest_bl_tdir}/CaseDocs/${tfile}"
    #fi
  done
  add_meshfile_to_docn_in "${tfile}" "${dest_bl_tdir}"
  # Create med_modelio.nml
  write_med_modelio "${dest_bl_tdir}/CaseDocs/med_modelio.nml"
  # Correct drv_in
  write_drv_in "${dest_bl_tdir}/CaseDocs/drv_in"
  # Update drv_in
  cat "${dest_bl_tdir}/CaseDocs/drv_in" | tail -221 | head -36 > "${dest_bl_tdir}/CaseDocs/temp"
  mv "${dest_bl_tdir}/CaseDocs/temp" "${dest_bl_tdir}/CaseDocs/drv_in"
  # Update docn_in
  if [ -f "${dest_bl_tdir}/CaseDocs/docn_in" ]; then
    cat "${dest_bl_tdir}/CaseDocs/docn_in" | head -7 > "${dest_bl_tdir}/CaseDocs/temp"
    mv "${dest_bl_tdir}/CaseDocs/temp" "${dest_bl_tdir}/CaseDocs/docn_in"
  fi
  # Create nuopc.runconfig (copy from store)
  if [ -f "${nuopc_store}/nuopc.runconfig.${nuopc_dir}" ]; then
    nfile="${nuopc_store}/nuopc.runconfig.${nuopc_dir}"
    cp "${nfile}" "${dest_bl_tdir}/CaseDocs/nuopc.runconfig"
  fi
  # Create nuopc.runseq
  if [ -f "${nuopc_store}/nuopc.runseq.${nuopc_dir}" ]; then
    nfile="${nuopc_store}/nuopc.runseq.${nuopc_dir}"
    cp "${nfile}" "${dest_bl_tdir}/CaseDocs/nuopc.runseq"
  fi
  # modify some <component>_modelio.nml files
  for component in atm glc ice lnd ocn rof wav; do
    fname="${dest_bl_tdir}/CaseDocs/${component}_modelio.nml"
    if [ -f "${fname}" ]; then
      sed -i -e 's/"pnetcdf/"netcdf/' ${fname}
    else
      perr "No ${component}_modelio.nml at '${fname}'"
    fi
  done
}

process_test() {
  # Args: <orig_bl> <driver_bl> <test dir> <driver>
  local orig_bl="${1}"
  local driver_bl="${2}"
  local test_dir="${3}"
  local driver="${4}"
  # Create baseline directory
  if [ ! -d "${driver_bl}" ]; then
    mkdir "${driver_bl}"
  fi
  # Process the baseline test directory
  clone_bl_test_dir "${orig_bl}" "${driver_bl}" "${test_dir}" "${driver}"
}

create_driver_baselines() {
  # Args: cam tag, baseline root, driver, compiler
  local ctag="${1}"      # CAM tag
  local bl_root="${2}"   # Original baseline root
  local dest_root="${3}" # Destination baseline root
  local driver="${4}"    # Driver to add (mct or nuopc)
  local comp="${5}"      # Compiler (nag, gnu, intel)
  local orig_bl
  local driver_bl
  local tdir
  local fmod
  if [ $# -ne 5 ]; then
    perr "Usage: ${0} <cam tag> <bl root> <dest root> <driver> <compiler>" 4
  fi
  orig_bl="${bl_root}/${ctag}_${comp}"
  driver_bl="${dest_root}/${ctag}_${driver}_${comp}"
  if  [ ! -d "${orig_bl}" ]; then
    orig_bl="${bl_root}/${ctag}"
    driver_bl="${dest_root}/${ctag}_${driver}"
  fi
  if  [ ! -d "${orig_bl}" ]; then
    perr "ERROR: Original baseline, \"${orig_bl}[_${comp}]\", must exist"
  fi
  if [ -d "${driver_bl}" ]; then
    echo "Deleting existing driver baselines"
    rm -r ${driver_bl}
  fi
  mkdir ${driver_bl}
  for tdir in `cd ${orig_bl}; ls`; do
    process_test "${orig_bl}" "${driver_bl}" "${tdir}" "${driver}"
  done
#  # modify esp_modelio.nml
#  fmod="s/pio_netcdf_format\ =\ [\"][\"]/pio_netcdf_format = \"64bit_offset\"/"
#  find ${driver_bl} -name esp_modelio.nml -exec sed -i -e "${fmod}" {} \;
  # modify drv_in
  fmod="s/case_name\([^.]*\)./case_name\1_V${driver}./"
  find ${driver_bl} -name drv_in -exec sed -i -e ${fmod} {} \;
  # Make softlink copies for any test in this baseline which keeps MCT support
  for bl_dir in ${MCT_CP}; do
    if [ -d "${orig_bl}/${bl_dir}" ]; then
      softlink_file "${bl_dir}" "${orig_bl}" "${driver_bl}"
    fi
  done
}

# Single argument is the tag to clone
# Optional second argument is an alternative destination for the new baselines
if [ $# -lt 1 -o $# -gt 2 ]; then
  perr "Usage: ${0} <cam tag> [ <dest root> ]"
fi

cam_tag="${1}"
if [ $# -gt 1 ]; then
  dest_root="${2}"
else
  dest_root=""
fi

if [ "${HOST:0:5}" == "izumi" ]; then
  if [ -z "${dest_root}" ]; then
    dest_root="${BL_ROOT}"
  fi
  # First, process nag
  create_driver_baselines ${cam_tag} ${BL_ROOT} ${dest_root} nuopc nag
  # Now, process gnu
  create_driver_baselines ${cam_tag} ${BL_ROOT} ${dest_root} nuopc gnu
elif [ "${HOST:0:8}" == "cheyenne" ]; then
  if [ -z "${dest_root}" ]; then
    dest_root="${BL_ROOT}"
  fi
  # Process intel tests
  create_driver_baselines ${cam_tag} ${BL_ROOT} ${dest_root} nuopc intel
else
  perr "ERROR: Unsupported host, \"${HOST}\"." 3
fi
