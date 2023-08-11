export interface Plate {
  barcode: string
}

export interface Well {
  location: string
  plate: Plate
}

export interface Sample {
  sample_number: string
  well: Well
  pseudoanonymized_id: string
}

export interface Strain {
  name: string
}

export interface Panel {
  name: string
  strain: Strain
}

export interface Substrain {
  name: string
  strain: Strain
}

export interface SampleCount {
  plate: Plate
  substrain: Substrain
  aligned: number
  length: number
  rpkm: number
  rpkm_proportions: number
  normcounts: number
  outlier: boolean
}

export interface DataItem {
  panel: string
  value: number | string | null
  strain: string
}

export interface Metadata {
  plate: Plate
  sample: Sample
  well: Well
  prescriber: string
  order_date: string
  ent_date: string
  treatment_type: string
  data: Array<DataItem>
}
