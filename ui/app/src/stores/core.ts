import {defineStore} from 'pinia'
import {api} from 'src/boot/axios'
import {ref} from 'vue'
import {Metadata, Plate, SampleCount, Substrain} from 'src/models/core'

export const useCoreStore = defineStore('core', () => {
  const sampleCounts = ref<SampleCount[]>([])
  const tableData = ref<SampleCount[]>([])
  const aggregate = ref<boolean>(false)
  const selected_barcode = ref<string | null>(null)
  const metadata = ref<Metadata[]>([])
  const plates = ref<Plate[]>([])
  const substrains = ref<Substrain[]>([])

  const toggleAggregate = () => {
    const mappedData = new Map()
    if (aggregate.value) {
      sampleCounts.value.forEach(item => {
        const strain = item.substrain.strain.name
        const substrain = item.substrain
        const existingData = mappedData.get(strain) || {
          aligned: 0,
          length: 0,
          rpkm: 0,
          rpkm_proportions: 0,
          normcounts: 0,
          outlier: false,
          plate: null,
          substrain: null,
        }

        const newData = {
          ...existingData,
          aligned: existingData.aligned + item.aligned,
          length: existingData.length + item.length,
          rpkm: existingData.rpkm + item.rpkm,
          rpkm_proportions: existingData.rpkm_proportions + item.rpkm_proportions,
          normcounts: existingData.normcounts + item.normcounts,
          outlier: item.outlier, // Update if needed
          plate: item.plate.barcode,
          substrain: substrain,
        }

        mappedData.set(strain, newData)
      })
      tableData.value = Array.from(mappedData.entries()).map(([key, value]) => ({
        strain: key,
        ...value,
      }))
    } else {
      tableData.value = sampleCounts.value
    }
  }

  const getPlates = async () => {
    try {
      const res = await api.get('/api/plates/')
      plates.value = res.data
    } catch (error) {
      console.error(error)
    }
  }

  const getSubstrains = async () => {
    try {
      const res = await api.get('/api/substrains/')
      substrains.value = res.data
    } catch (error) {
      console.error(error)
    }
  }

  const getPlateData = async (barcode: string | null = null, substrain: string | null = null) => {
    let baseUrl = '/api/samplecounts/'

    if (barcode) {
      baseUrl += `?plate__barcode=${barcode}`
      selected_barcode.value = barcode

      if (substrain) {
        baseUrl += `&substrain__name=${substrain}`
      }

      try {
        const res1 = await api.get(baseUrl)
        sampleCounts.value = res1.data
        tableData.value = res1.data
        const res2 = await api.get(`/api/metadata/?plate__barcode=${barcode}`)
        metadata.value = res2.data
      } catch (error) {
        console.error(error)
      }
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
    toggleAggregate,
  }
})
