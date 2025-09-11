import {defineStore} from 'pinia'
import {api} from 'src/boot/axios'
import {ref} from 'vue'
import {Metadata, Plate, SampleCount, Substrain, Sample} from 'src/models/core'

export const useCoreStore = defineStore('core', () => {
  const sampleCounts = ref<SampleCount[]>([])
  const tableData = ref<SampleCount[]>([])
  const aggregate = ref<boolean>(false)
  const selected_barcode = ref<string | null>(null)
  const metadata = ref<Metadata[]>([])
  const plates = ref<Plate[]>([])
  const substrains = ref<Substrain[]>([])
  const samples = ref<Sample[]>([])
  const selected_plate_barcode = ref<string | null>(null)
  const selected_plate = ref<Plate | null>(null)
  const selected_substrain = ref<string | null>(null)
  const selected_sample_id = ref<string | null>(null)
  const selected_sample = ref<Sample | null>(null)

// Helper: make DRF `next` safe for the axios `baseURL`
const toRelative = (u: string) => u.replace(/^https?:\/\/[^/]+/, '');

// Helper: fetch all pages (DRF-style or plain arrays)
const fetchAllPages = async (initialUrl: string) => {
  const all: any[] = [];
  let url: string | null = initialUrl;
  while (url) {
    const res: any = await api.get(url);
    const data: any = res.data;
    if (Array.isArray(data)) {   // Non-paginated list endpoint
      all.push(...data);
      url = null;
    } else if (data && Array.isArray(data.results)) { // DRF paginated
      all.push(...data.results);
      url = data.next ? toRelative(data.next) : null;
    } else {
      return data; // Non-list object; return as-is
    }
  }
  return all;
};

  const getSelectedSample = async () => {
    try {
      const res = await api.get(`/api/samples/?pseudonymized_id=${selected_sample_id.value}`)
      selected_sample.value = res.data
    } catch (error) {
      console.error(error)
    }
  }

  const filterCountDataBySample = async (pseudonymized_id = '') => {
    if (pseudonymized_id !== '') {
      tableData.value = sampleCounts.value.filter(
        item => item.sample.pseudonymized_id === pseudonymized_id
      )
    } else {
      tableData.value = sampleCounts.value
    }

    await getSelectedSample()
  }

  const aggregateData = () => {
    aggregate.value = true
    const mappedData = new Map()
      tableData.value.forEach(item => {
        const strain_sampleId = item.substrain.strain.name + '__' + item.sample.pseudonymized_id
        // const strain = item.substrain.strain.name
        const substrain = item.substrain
        const existingData = mappedData.get(strain_sampleId) || {
          aligned: 0,
          length: 0,
          rpkm: 0,
          rpkm_proportions: 0,
          normcounts: 0,
          outlier: false,
          coverage_threshold: 0,
          coverage: 0,
          coverage_status: '',
          readnum_threshold: 0,
          readnum_status: '',
          percentile_threshold: '',
          plate: null,
          substrain: null,
          pseudonymized_id: null,
        }

        const newData = {
          ...existingData,
          aligned: existingData.aligned + item.aligned,
          length: existingData.length + item.length,
          rpkm: existingData.rpkm + item.rpkm,
          rpkm_proportions: existingData.rpkm_proportions + item.rpkm_proportions,
          normcounts: existingData.normcounts + item.normcounts,
          coverage_threshold:  item.coverage_threshold,
          coverage: existingData.coverage + item.coverage,
          readnum_threshold:  item.readnum_threshold,
          percentile_threshold: item.percentile_threshold,
          outlier: item.outlier,
          plate: item.plate.barcode,
          substrain: substrain,
          pseudonymized_id: item.sample.pseudonymized_id,
        }

        mappedData.set(strain_sampleId, newData)
      })
    console.log(mappedData)

      tableData.value = Array.from(mappedData.entries()).map(([key, value]) => ({
        strain: key.split('__')[0],
        pseudonymized_id: key.split('__')[1],
        ...value,
      }))

  }

  const cancelAggregate = () => {
    aggregate.value = false
    tableData.value = sampleCounts.value
  }

const getPlates = async (barcode: string | null = null) => {
  try {
    if (barcode) { // Fetch a specific plate by barcode
      const res = await api.get(`/api/plates/?barcode=${barcode}`);
      selected_plate.value = res.data[0] || null;
    } else { // Fetch all plates using the helper
      plates.value = await fetchAllPages('/api/plates/');
    }
  } catch (error) {
    console.error('Error fetching plates:', error);
  }
};

const getSubstrains = async () => {
  try {
    substrains.value = await fetchAllPages('/api/substrains/');
  } catch (error) {
    console.error('Error fetching substrains:', error);
  }
};

const getSamplesByPlate = async (barcode: string) => {
  try {
    samples.value = await fetchAllPages(`/api/samples/?plate__barcode=${barcode}`);
  } catch (error) {
    console.error('Error fetching samples by plate:', error);
  }
};

const getPlateData = async (barcode: string | null = null, substrain: string | null = null) => {
  try {
    if (!barcode) return;
    selected_barcode.value = barcode;
    await getPlates(barcode);
    const params = new URLSearchParams({ 'plate__barcode': barcode }); // build query parameters with URLSearchParams to handle encoding
    if (substrain) {
      params.append('substrain__name', substrain);
    }
    const [counts, meta] = await Promise.all([ // fetch all sample counts and metadata in parallel
      fetchAllPages(`/api/samplecounts/?${params.toString()}`),
      fetchAllPages(`/api/metadata/?plate__barcode=${barcode}`)
    ]);
    sampleCounts.value = counts;
    tableData.value = counts;
    metadata.value = meta;

  } catch (error) {
    console.error('Error fetching plate data:', error);
  }
};

  const downloadFile = async (path: string) => {
    try {
      const encodedPath = encodeURIComponent(path)
      const res = await api.get(`/api/download/${encodedPath}/`, {
        responseType: 'blob',
      })
      const url = window.URL.createObjectURL(new Blob([res.data]))
      const link = document.createElement('a')
      link.href = url
      link.setAttribute('download', path.split('/').pop() || '')
      document.body.appendChild(link)
      link.click()
    } catch (error) {
      console.error('Error downloading the file:', error)
    }
  }

  return {
    getPlateData,
    sampleCounts,
    selected_barcode,
    metadata,
    getPlates,
    plates,
    getSubstrains,
    substrains,
    tableData,
    aggregate,
    aggregateData,
    getSamplesByPlate,
    samples,
    filterCountDataBySample,
    selected_plate_barcode,
    selected_substrain,
    selected_sample,
    selected_plate,
    downloadFile,
    selected_sample_id,
    cancelAggregate,
  }
})
