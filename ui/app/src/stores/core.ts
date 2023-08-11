import {defineStore} from 'pinia'
import {api} from 'src/boot/axios'
import {ref} from 'vue'
import {SampleCount, Metadata, Plate, Substrain} from 'src/models/core'

export const useCoreStore = defineStore('core', () => {
  const sampleCounts = ref<SampleCount[]>([])
  const selected_barcode = ref<string | null>(null)
  const metadata = ref<Metadata[]>([])
  const plates = ref<Plate[]>([])
  const substrains = ref<Substrain[]>([])

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
      baseUrl += `?barcode=${barcode}`
      selected_barcode.value = barcode

      if (substrain) {
        baseUrl += `&substrain=${substrain}`
      }

      try {
        const res1 = await api.get(baseUrl)
        sampleCounts.value = res1.data
        const res2 = await api.get(`/api/metadata/?barcode=${barcode}`)
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
  }
})
